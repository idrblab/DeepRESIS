from __future__ import annotations

import pytest

from deepresis.cli import main


def test_cli_help(capsys):
    assert main([]) == 0
    out = capsys.readouterr().out
    assert "deepresis" in out


def test_cli_predict_help(capsys):
    with pytest.raises(SystemExit) as exc_info:
        main(["predict", "--help"])
    assert exc_info.value.code == 0
    out = capsys.readouterr().out
    assert "--fasta" in out
