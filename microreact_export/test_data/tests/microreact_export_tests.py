import argparse
import base64
import importlib.util
import json
import logging
import os
from pathlib import Path
from uuid import UUID

import pytest


@pytest.fixture(scope="session")
def micro_react():
    env_path = os.getenv("MICROREACT_EXPORT_PATH")
    if env_path:
        module_path = Path(env_path).expanduser().resolve()
        if module_path.is_file() == False:
            raise FileNotFoundError(f"MICROREACT_EXPORT_PATH does not exist: {module_path}")

    spec = importlib.util.spec_from_file_location("microreact_export_module", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def logger():
    logger = logging.getLogger("microreact_export_tests")
    logger.setLevel(logging.DEBUG)
    if not logger.handlers:
        logger.addHandler(logging.NullHandler())
    return logger


class FakeResponse:
    def __init__(self, status_code=200, payload=None, text="OK", raise_http=False):
        self.status_code = status_code
        self._payload = payload or {}
        self.text = text
        self._raise_http = raise_http

    def raise_for_status(self):
        if self._raise_http:
            import requests

            raise requests.exceptions.HTTPError(self.text)

    def json(self):
        return self._payload


def decode_blob(blob: str) -> str:
    encoded = blob.split(",", 1)[1]
    return base64.b64decode(encoded).decode("utf-8")


def test_logger(micro_react, tmp_path):
    log_file = tmp_path / "test.log"
    logger = micro_react.setup_logger(str(log_file), verbose=True)

    logger.debug("debug-message")
    for handler in logger.handlers:
        try:
            handler.flush()
        except Exception:
            pass

    assert log_file.exists()
    assert "debug-message" in log_file.read_text(encoding="utf-8")
    assert logger.level == logging.DEBUG


def test_parse_metadata_tsv_basic(micro_react, tmp_path, logger):
    tsv = tmp_path / "meta.tsv"
    tsv.write_text("sample_id\tcountry\ns1\tUSA\ns2\tCAN\n", encoding="utf-8")

    csv_data, headers, id_col, dates_validated = micro_react.parse_metadata_tsv(
        tsv=tsv,
        id_column="sample_id",
        logger=logger,
        remove_file_columns=False,
        selected_columns=None,
        date_column=None,
    )

    assert "sample_id,country" in csv_data
    assert headers == ["sample_id", "country"]
    assert id_col == "sample_id"
    assert dates_validated is False


def test_parse_metadata_tsv_remove_file_columns_and_date(micro_react, tmp_path, logger):
    tsv = tmp_path / "meta.tsv"
    tsv.write_text(
        "sample_id\turl\tcollection_date\n"
        "S1\thttps://x/y\t2023-01-01\n"
        "S2\tftp://a/b\tnot-a-date\n",
        encoding="utf-8",
    )

    csv_data, headers, id_col, dates_validated = micro_react.parse_metadata_tsv(
        tsv=tsv,
        id_column="sample_id",
        logger=logger,
        remove_file_columns=True,
        selected_columns=None,
        date_column="collection_date",
    )

    assert "url" not in headers
    assert "collection_date_validated" in headers
    assert id_col == "sample_id"
    assert dates_validated is True
    assert "sample_id,collection_date,collection_date_validated" in csv_data


def test_parse_metadata_tsv_defaults_id_column_if_missing(micro_react, tmp_path, logger):
    tsv = tmp_path / "meta.tsv"
    tsv.write_text("first\tsecond\nA\tB\n", encoding="utf-8")

    _, headers, id_col, _ = micro_react.parse_metadata_tsv(
        tsv=tsv,
        id_column="not_present",
        logger=logger,
    )

    assert headers[0] == "first"
    assert id_col == "first"


def test_encode_returns_prefixed_blob_and_uuid(micro_react, logger):
    encoded, file_id = micro_react.encode("hello", logger)
    assert encoded.startswith("data:application/octet-stream;base64,")
    UUID(file_id)


def test_add_entry_new_and_update(micro_react, logger):
    project = {}
    project = micro_react.add_entry(project, "files", {"a": 1}, logger)
    assert project["files"] == {"a": 1}

    project = micro_react.add_entry(project, "files", {"b": 2}, logger)
    assert project["files"] == {"a": 1, "b": 2}


def test_create_tree_entry_valid(micro_react, tmp_path, logger):
    tree = tmp_path / "t1.nwk"
    tree.write_text("(A:1,B:1);", encoding="utf-8")

    files_dict, tree_dict, tree_names = micro_react.create_tree_entry(
        tree_files=[tree], id_column="sample_id", set_id="setX", logger=logger
    )

    assert len(files_dict) == 1
    assert len(tree_dict) == 1
    assert len(tree_names) == 1

    name = "setX_t1"
    tree_id = tree_names[name]
    assert tree_dict[name]["labelField"] == "sample_id"
    assert files_dict[tree_id]["type"] == "tree"


def test_create_tree_entry_unsupported_type(micro_react, tmp_path, logger):
    bad = tmp_path / "t1.txt"
    bad.write_text("x", encoding="utf-8")

    with pytest.raises(RuntimeError, match="No valid tree files created"):
        micro_react.create_tree_entry([bad], "id", "setX", logger)


def test_create_metadata_entry_with_and_without_timeline(micro_react, logger):
    entry, timeline, entry_id = micro_react.create_metadata_entry(
        metadata_parsed="a,b\n1,2\n",
        dates_validated=True,
        logger=logger,
        date_column="date",
    )
    assert entry_id in entry
    assert "date" in timeline
    assert timeline["date"]["valueField"] == "date_validated"

    _, timeline2, _ = micro_react.create_metadata_entry(
        metadata_parsed="a,b\n1,2\n",
        dates_validated=False,
        logger=logger,
        date_column="date",
    )
    assert timeline2 == {}


def test_create_map_entry(micro_react, logger):
    map_entry = micro_react.create_map_entry(logger)
    assert len(map_entry) == 1
    map_id = next(iter(map_entry))
    assert map_entry[map_id]["type"] == "mapbox"
    assert map_entry[map_id]["latitudeField"] == "latitude"


def test_create_matrix_entry(micro_react, tmp_path, logger):
    matrix = tmp_path / "dist.csv"
    matrix.write_text("snp-dists 0.8.2,A,B\nA,0,1\nB,1,0\n", encoding="utf-8")

    matrix_file_entry, matrix_entry, matrix_name = micro_react.create_matrix_entry(
        matrix_files=[matrix],
        set_id="setM",
        logger=logger,
    )

    assert matrix_name == "setM_dist"
    assert len(matrix_file_entry) == 1
    assert len(matrix_entry) == 1

    matrix_id = next(iter(matrix_file_entry))
    blob = matrix_file_entry[matrix_id]["blob"]
    decoded_csv = decode_blob(blob)
    assert decoded_csv.splitlines()[0].startswith(",")


def test_submit_microreact_project_success(micro_react, logger, monkeypatch):
    captured = {}

    def fake_post(url, headers=None, json=None):
        captured["url"] = url
        captured["headers"] = headers
        captured["json"] = json
        return FakeResponse(status_code=200, payload={"id": "proj123"})

    monkeypatch.setattr(micro_react.requests, "post", fake_post)

    result = micro_react.submit_microreact_project(
        project_input={"x": 1},
        access_token="abc",
        logger=logger,
        restricted_access=True,
    )

    assert result["id"] == "proj123"
    assert captured["url"].endswith("/create?access=private")
    assert captured["headers"]["Access-Token"] == "abc"


def test_submit_microreact_project_http_error_raises(micro_react, logger, monkeypatch):
    def fake_post(url, headers=None, json=None):
        return FakeResponse(status_code=500, payload={}, text="boom", raise_http=True)

    monkeypatch.setattr(micro_react.requests, "post", fake_post)

    with pytest.raises(Exception):
        micro_react.submit_microreact_project(
            project_input={"x": 1},
            access_token="abc",
            logger=logger,
            restricted_access=False,
        )


@pytest.mark.xfail(reason="Current function returns undefined variable when access_token is empty, end to end functionality" \
    "provides that this cannot occur. Microreact submission won't occur unless the token is provided." \
    "This function can be made more robust with extra error handling, however, and this case will be covered.")
def test_submit_microreact_project_without_token(micro_react, logger):
    micro_react.submit_microreact_project(project_input={"x": 1}, access_token="", logger=logger)


def test_update_microreact_project_minimal(micro_react, logger, monkeypatch):
    existing = {
        "datasets": {},
        "files": {},
        "tables": {},
        "trees": {},
        "matrices": {},
        "maps": {},
        "timelines": {},
    }

    sent = {}

    def fake_get(url, headers=None):
        return FakeResponse(status_code=200, payload=existing)

    def fake_post(url, headers=None, json=None):
        sent["url"] = url
        sent["json"] = json
        return FakeResponse(status_code=200, payload={"ok": True})

    monkeypatch.setattr(micro_react.requests, "get", fake_get)
    monkeypatch.setattr(micro_react.requests, "post", fake_post)

    resp, payload_str = micro_react.update_microreact_project(
        project_url="p123",
        access_token="tok",
        id_column="sample_id",
        set_id="set1",
        metadata=None,
        date_column=None,
        remove_file_columns=True,
        tree_files=None,
        matrix_files=None,
        logger=logger,
    )

    assert resp.status_code == 200
    payload = json.loads(payload_str)
    assert payload == sent["json"]


def test_update_microreact_project_with_metadata_tree_and_matrix(micro_react, logger, tmp_path, monkeypatch):
    metadata = tmp_path / "meta.tsv"
    metadata.write_text(
        "sample_id\tlatitude\tlongitude\tdate\nS1\t1.0\t2.0\t2024-01-01\n",
        encoding="utf-8",
    )
    tree = tmp_path / "phylo.nwk"
    tree.write_text("(S1:1,S2:1);", encoding="utf-8")
    matrix = tmp_path / "dist.csv"
    matrix.write_text("x,S1\nS1,0\n", encoding="utf-8")

    tree_name = "set1_phylo"
    existing = {
        "datasets": {},
        "files": {"old_tree_file_id": {"id": "old_tree_file_id"}},
        "tables": {},
        "trees": {tree_name: {"file": "old_tree_file_id"}},
        "matrices": {},
        "maps": {},
        "timelines": {},
    }

    sent = {}

    def fake_get(url, headers=None):
        return FakeResponse(status_code=200, payload=existing)

    def fake_post(url, headers=None, json=None):
        sent["json"] = json
        return FakeResponse(status_code=200, payload={"ok": True})

    monkeypatch.setattr(micro_react.requests, "get", fake_get)
    monkeypatch.setattr(micro_react.requests, "post", fake_post)

    micro_react.update_microreact_project(
        project_url="p123",
        access_token="tok",
        id_column="sample_id",
        set_id="set1",
        metadata=metadata,
        date_column="date",
        remove_file_columns=True,
        tree_files=[tree],
        matrix_files=[matrix],
        logger=logger,
    )

    payload = sent["json"]
    assert payload["maps"] is not None
    assert len(payload["datasets"]) == 1
    assert len(payload["tables"]) == 1
    assert "date" in payload["timelines"]
    assert payload["trees"][tree_name]["labelField"] == "sample_id"
    assert payload["matrices"]


def test_create_microreact_project_full(micro_react, logger, tmp_path):
    metadata = tmp_path / "meta.tsv"
    metadata.write_text(
        "sample_id\tlatitude\tlongitude\tdate\tcountry\n"
        "S1\t1.1\t2.2\t2024-01-01\tUSA\n",
        encoding="utf-8",
    )
    tree = tmp_path / "tree.nwk"
    tree.write_text("(S1:1,S2:1);", encoding="utf-8")
    matrix = tmp_path / "dist.csv"
    matrix.write_text("x,S1\nS1,0\n", encoding="utf-8")

    project_json = micro_react.create_microreact_project(
        metadata=metadata,
        id_column="sample_id",
        set_id="set9",
        matrix_files=[matrix],
        date_column="date",
        tree_files=[tree],
        logger=logger,
        remove_file_columns=True,
        project_name="My Project",
        selected_columns=["sample_id", "date", "country", "missing_col"],
    )

    payload = json.loads(project_json)
    assert payload["meta"]["name"] == "My Project"
    assert payload["datasets"]
    assert payload["files"]
    assert payload["tables"]
    assert payload["trees"]
    assert payload["matrices"]
    assert payload["maps"]
    assert payload["timelines"]


def test_create_microreact_project_requires_metadata(micro_react, logger):
    with pytest.raises(ValueError, match="Metadata file path must be provided"):
        micro_react.create_microreact_project(
            metadata=None,
            id_column="sample_id",
            set_id="set9",
            matrix_files=None,
            date_column=None,
            tree_files=None,
            logger=logger,
        )


def test_main_create_path_writes_project_file(micro_react, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def fake_parse_args(self):
        return argparse.Namespace(
            project_name="projA",
            project_url=None,
            set_id="set1",
            metadata_tsv=Path("meta.tsv"),
            matrix_files=None,
            tree_files=None,
            selected_columns=None,
            access_token=None,
            restricted_access=False,
            remove_file_columns=False,
            id_column="sample_id",
            date_column=None,
            log="run.log",
            verbose=False,
        )

    monkeypatch.setattr(micro_react.argparse.ArgumentParser, "parse_args", fake_parse_args)
    monkeypatch.setattr(micro_react, "setup_logger", lambda *a, **k: logging.getLogger("main_test"))
    monkeypatch.setattr(micro_react, "create_microreact_project", lambda **kwargs: '{"k":"v"}')

    micro_react.main()

    out_file = tmp_path / "projA_input.microreact"
    assert out_file.exists()
    assert json.loads(out_file.read_text(encoding="utf-8")) == {"k": "v"}


def test_main_update_path_writes_project_file(micro_react, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def fake_parse_args(self):
        return argparse.Namespace(
            project_name="projB",
            project_url="https://microreact.org/project/x",
            set_id="set1",
            metadata_tsv=None,
            matrix_files=None,
            tree_files=None,
            selected_columns=None,
            access_token="tok",
            restricted_access=False,
            remove_file_columns=False,
            id_column="sample_id",
            date_column=None,
            log="run.log",
            verbose=False,
        )

    monkeypatch.setattr(micro_react.argparse.ArgumentParser, "parse_args", fake_parse_args)
    monkeypatch.setattr(micro_react, "setup_logger", lambda *a, **k: logging.getLogger("main_test"))
    monkeypatch.setattr(micro_react, "update_microreact_project", lambda **kwargs: (FakeResponse(), '{"updated":true}'))

    micro_react.main()

    out_file = tmp_path / "projB_input.microreact"
    assert out_file.exists()
    assert json.loads(out_file.read_text(encoding="utf-8")) == {"updated": True}