import logging
import sys
from pathlib import Path


def setup_logging(log_file: Path, mode: str = "w") -> None:
    """Logs to stdout and to log_file"""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.StreamHandler(sys.stdout), logging.FileHandler(log_file, mode=mode)],
    )
