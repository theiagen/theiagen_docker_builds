#!/usr/bin/env python3
"""Check that software installed in a container satisfies a requirements txt file.

Each line of the requirements file is <software><op><version>, op one of
== != >= <= > < ~=  (e.g. samtools>=1.17). Blank lines, # comments, and names
without an operator are ignored. Software not found by any manager is skipped.

Usage:
  check_dependencies.py deps.txt --image <image>          # throwaway container
  check_dependencies.py deps.txt --container <running>    # docker exec
  check_dependencies.py deps.txt                          # local environment
"""
import argparse
import operator
import re
import subprocess
import sys
from contextlib import contextmanager, nullcontext

# string operator to python operator dict
OPS = {"==": operator.eq, "!=": operator.ne, ">=": operator.ge,
       "<=": operator.le, ">": operator.gt, "<": operator.lt}
# regex for splitting software, operator, version
REQ_RE = re.compile(r"^\s*([A-Za-z0-9_.+-]+?)\s*(~=|==|!=|>=|<=|>|<)\s*([^\s,;]+)")
# multiple executables to check conda package manager
CONDA_EXES = ("conda", "mamba", "micromamba", '"$MAMBA_EXE"')


def _norm(name):
    return name.strip().lower().replace("_", "-")


def _vkey(version):
    """Numeric tuple from a version, dropping epochs (1:) and suffixes (-1ubuntu, +dfsg)."""
    # convert version to tuple
    core = re.split(r"[-+~]", version.split(":", 1)[-1])[0]
    # enforce the versions to be tuples of integers
    # may cause edge-case failures where 1.2.3a1 is actually < 1.2.3, but reported as after
    # because the alpha is stripped
    return tuple(int(x) for x in re.findall(r"\d+", core))


def _parse(text, pattern):
    return {_norm(m.group(1)): m.group(2) for m in re.finditer(pattern, text, re.M)}


def version_satisfies(found, op, required):
    """Query and check the versions; pads versions with 0s for comparison
    and handles suffixes + epochs"""
    a, b = _vkey(found), _vkey(required)
    k, n = len(b), max(len(a), len(b))
    # pad the version with 0s so they are directly comparable
    a, b = a + (0,) * (n - len(a)), b + (0,) * (n - k)
    if op == "~=":
        return a >= b and a[:k - 1] == b[:k - 1]
    return OPS[op](a, b)


def parse_requirements(path):
    """Return {software: (operator, version)} from a requirements txt file."""
    reqs = {}
    with open(path) as fh:
        for line in fh:
            # extract the software, operator, and version
            m = REQ_RE.match(line.split("#", 1)[0])
            if m:
                reqs[_norm(m.group(1))] = (m.group(2), m.group(3))
    return reqs


# --- command runners: run(cmd) -> stdout ("" on any failure) ---------------
def make_runner(prefix=()):
    """Run `sh -c cmd` behind an optional prefix (e.g. ["docker", "exec", id])."""
    def run(cmd):
        try:
            return subprocess.run(list(prefix) + ["sh", "-c", cmd], stdout=subprocess.PIPE,
                                  stderr=subprocess.DEVNULL, universal_newlines=True,
                                  timeout=600).stdout
        # do not fail if the software fails to execute, i.e. if the package manager does not
        # exist
        except (OSError, subprocess.SubprocessError):
            return ""
    return run


def exec_runner(container, engine="docker"):
    """Runner for an already-running container; raises if it is not reachable."""
    prefix = [engine, "exec", container]
    subprocess.run(prefix + ["sh", "-c", "true"], check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return make_runner(prefix)


@contextmanager
def image_runner(image, engine="docker"):
    """Start a throwaway container from `image`, yield a runner for it, then remove it."""
    cid = subprocess.check_output([engine, "run", "-d", "--rm", "--entrypoint", "tail",
                                   image, "-f", "/dev/null"], universal_newlines=True).strip()
    try:
        yield make_runner([engine, "exec", cid])
    finally:
        subprocess.run([engine, "rm", "-f", cid],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


# --- package managers: query(run) -> {software: version} -------------------
def query_apt(run):
    return _parse(run(r"dpkg-query -W -f='${Package}\t${Version}\n' 2>/dev/null"),
                  r"^(\S+)\t(\S+)$")


def query_pip(run):
    cmds = ("pip3", "pip", "python3 -m pip", "python -m pip")
    out = run("; ".join(c + " list --format=freeze 2>/dev/null" for c in cmds))
    return _parse(out, r"^([A-Za-z0-9_.-]+)==(\S+)$")


def query_micromamba(run):
    """Packages in every environment visible to conda, mamba, and micromamba."""
    found = {}
    for exe in CONDA_EXES:
        for prefix in re.findall(r"(/\S*)\s*$", run(exe + " env list 2>/dev/null"), re.M):
            out = run("%s list -p '%s' 2>/dev/null" % (exe, prefix))
            found.update(_parse(out, r"^\s*([A-Za-z0-9_.+-]+)\s+(\d\S*)\s"))
    return found


MANAGERS = {"apt": query_apt, "pip": query_pip, "micromamba": query_micromamba}


def check_dependencies(requirements, run=None, managers=None):
    """Enforce requirements (path or parse_requirements() dict) against each manager.

    Returns {software: {manager: version}} for software that was found; raises
    ValueError listing every installed version that violates its constraint.
    """
    if isinstance(requirements, str):
        requirements = parse_requirements(requirements)
    run, managers = run or make_runner(), managers or MANAGERS
    found, failures = {}, []
    for manager, query in managers.items():
        installed = query(run)
        for name, (op, required) in requirements.items():
            if name not in installed:
                continue
            version = found.setdefault(name, {})[manager] = installed[name]
            if not version_satisfies(version, op, required):
                failures.append("%s %s (%s) does not satisfy %s%s"
                                % (name, version, manager, op, required))
    if failures:
        raise ValueError("Unmet dependency requirements:\n  " + "\n  ".join(failures))
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("requirements", help="txt file of <software><op><version> lines")
    target = ap.add_mutually_exclusive_group()
    target.add_argument("-i", "--image", help="image to start a throwaway container from")
    target.add_argument("-c", "--container", help="name/ID of a running container")
    ap.add_argument("-e", "--engine", default="docker", help="container engine [docker]")
    ap.add_argument("-m", "--managers", nargs="+", choices=list(MANAGERS), default=list(MANAGERS))
    args = ap.parse_args()

    reqs = parse_requirements(args.requirements)
    managers = {m: MANAGERS[m] for m in args.managers}
    if args.image:
        ctx = image_runner(args.image, args.engine)
    else:
        ctx = nullcontext(exec_runner(args.container, args.engine)
                          if args.container else make_runner())
    try:
        with ctx as run:
            found = check_dependencies(reqs, run, managers)
    except ValueError as e:
        sys.exit(str(e))
    for name, (op, required) in sorted(reqs.items()):
        hits = ", ".join("%s %s" % (v, m) for m, v in found.get(name, {}).items())
        print("%-8s %s%s%s  %s" % ("OK" if hits else "MISSING", name, op, required, hits))


if __name__ == "__main__":
    main()
