"""Test functions of util.py ."""

from cobra.util import AutoVivification, format_long_string, show_versions


def test_format_long_string_with_long_string() -> None:
    """Test functionality by formatting a long string."""
    long_string = (
        "This is a really long string, but is it long enough."
        "I hope it is long enough so that format_long_string() function works."
    )
    expected_output = "This is a really long string, but is it long en..."
    assert expected_output == format_long_string(long_string)


def test_format_long_string_with_short_string() -> None:
    """Test functionality by formatting a short string."""
    short_string = "This is short string."
    assert short_string == format_long_string(short_string)


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
    assert lines[1].startswith("System Information")
    assert lines[2].startswith("==================")
    assert lines[3].startswith("OS")
    assert lines[4].startswith("OS-release")
    assert lines[5].startswith("Python")

    assert lines[7].startswith("Package Versions")
    assert lines[8].startswith("================")
