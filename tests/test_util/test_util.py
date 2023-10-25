"""Test functions of util.py ."""

import pytest

from cobra.util import AutoVivification, format_long_string, show_versions


@pytest.mark.parametrize(
    "input_string, expected_string",
    [
        (
            (
                "This is a really long string, but is it long enough. "
                "I hope it is long enough so that format_long_string() function works."
            ),
            "This is a really long string, but is it long...",
        ),
        ("This is short string.", "This is short string."),
    ],
)
def test_format_long_string(input_string: str, expected_string: str) -> None:
    """Test functionality of format long string."""
    assert expected_string == format_long_string(input_string)


def test_autovivification() -> None:
    """Test proper functionality of autovivification."""
    test_data = AutoVivification()
    test_data["a"]["b"] = 1
    test_data["c"]["d"] = 2
    assert test_data["a"] == {"b": 1}
    assert test_data["c"] == {"d": 2}
    assert test_data["a"]["b"] == 1
    assert test_data["c"]["d"] == 2


def test_show_versions(capsys) -> None:
    """Test output of dependency information."""
    show_versions()
    captured = capsys.readouterr()
    lines = captured.out.split("\n")
    assert lines[1].startswith("Package Information")
    assert lines[2].startswith("------------------")
    assert lines[3].startswith("cobra")

    assert lines[5].startswith("Dependency Information")
    assert lines[6].startswith("------------------")
