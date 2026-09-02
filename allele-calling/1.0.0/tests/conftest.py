import os
from pathlib import Path

import pytest


@pytest.fixture
def change_test_dir_to_temp(tmp_path: Path, request: pytest.FixtureRequest):
    os.chdir(tmp_path)
    yield
    os.chdir(request.config.invocation_params.dir)
