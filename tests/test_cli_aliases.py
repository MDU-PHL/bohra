from click.testing import CliRunner

import bohra.bohra as bohra_cli
from bohra.commands.test import bohratest


def test_test_command_accepts_threads_alias(monkeypatch, tmp_path):
    captured = {}

    def fake_run_tests(cpus, shovill_ram, wdir):
        captured["cpus"] = cpus
        captured["shovill_ram"] = shovill_ram
        captured["wdir"] = wdir

    monkeypatch.setattr(bohratest, "run_tests", fake_run_tests)
    result = CliRunner().invoke(
        bohratest.test,
        ["--threads", "3", "--shovill_ram", "8", "--wdir", str(tmp_path)],
    )

    assert result.exit_code == 0
    assert captured == {"cpus": 3, "shovill_ram": 8, "wdir": str(tmp_path)}


def test_run_command_accepts_threads_alias(monkeypatch):
    captured = {}

    def fake_run_bohra(pipeline, kwargs):
        captured["pipeline"] = pipeline
        captured["kwargs"] = kwargs

    monkeypatch.setattr(bohra_cli, "run_bohra", fake_run_bohra)
    monkeypatch.setattr(
        bohra_cli.sys,
        "argv",
        ["bohra", "run", "basic", "--threads", "4", "--input_file", "input.tsv"],
    )
    result = CliRunner().invoke(
        bohra_cli.cli,
        ["run", "basic", "--threads", "4", "--input_file", "input.tsv"],
    )

    assert result.exit_code == 0
    assert captured["pipeline"] == "basic"
    assert captured["kwargs"]["cpus"] == 4
