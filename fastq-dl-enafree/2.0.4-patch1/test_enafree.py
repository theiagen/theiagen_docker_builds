"""Checks that the fastq-dl 2.0.4 ENA-free patch is in place. No network access is needed."""

from unittest import mock

import pandas as pd
import requests
from fastq_dl import fastq_dl


def test_ena_connection_error_returns_failure():
    connection_error = requests.exceptions.ConnectionError("reset")
    with mock.patch.object(requests, "get", side_effect=connection_error) as get:
        success, payload = fastq_dl.get_ena_metadata("SRX000001")
    assert success is False
    assert payload == [None, "reset"]
    assert get.call_args.kwargs["timeout"] == 30


def test_sra_metadata_skips_pysradb_ena_lookup():
    with mock.patch.object(fastq_dl, "SRAweb") as sraweb:
        sraweb.return_value.search_sra.return_value = None
        fastq_dl.get_sra_metadata("SRX000001")
        ena_lookup = sraweb.return_value.fetch_ena_fastq("SRP000001")
    assert isinstance(ena_lookup, pd.DataFrame)
    assert ena_lookup.empty


def test_provider_sra_does_not_query_ena():
    sra_records = [{"run_accession": "SRR000001"}]
    with (
        mock.patch.object(
            fastq_dl, "get_ena_metadata", side_effect=AssertionError("ENA queried")
        ),
        mock.patch.object(
            fastq_dl, "get_sra_metadata", return_value=[True, sra_records]
        ),
    ):
        provider, records = fastq_dl.get_run_info(
            "SRX000001", "SRX000001", "sra", False, max_attempts=1, sleep=0
        )
    assert provider == fastq_dl.SRA
    assert records == sra_records


def test_provider_sra_failure_does_not_fall_back_to_ena():
    with (
        mock.patch.object(
            fastq_dl, "get_ena_metadata", side_effect=AssertionError("ENA queried")
        ),
        mock.patch.object(fastq_dl, "get_sra_metadata", return_value=[False, []]),
    ):
        try:
            fastq_dl.get_run_info(
                "SRX000001", "SRX000001", "sra", False, max_attempts=1, sleep=0
            )
        except SystemExit as error:
            assert error.code == 1
        else:
            raise AssertionError("get_run_info did not exit after the SRA failure")


if __name__ == "__main__":
    for name, test in sorted(globals().items()):
        if name.startswith("test_"):
            test()
            print(f"PASS {name}")
