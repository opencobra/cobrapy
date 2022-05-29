from cobra.core.metadata.keyvaluepairs import KeyValueEntry, KeyValuePairs


EXAMPLE_URL = "https://tinyurl.com/ybyr7b62"


def test_keyvalueentry():
    keyvaluedict = KeyValueEntry.from_data(
        {
            "id": "KV_id",
            "name": "abc_xyz",
            "key": "keyX",
            "value": "45",
            "uri": EXAMPLE_URL,
        }
    )
    assert isinstance(keyvaluedict, KeyValueEntry)
    assert keyvaluedict.id == "KV_id"
    assert keyvaluedict.name == "abc_xyz"
    assert keyvaluedict.key == "keyX"
    assert keyvaluedict.value == "45"
    assert keyvaluedict.uri == EXAMPLE_URL


def test_keyvalueentry_empty_uri():
    keyvaluedict = KeyValueEntry.from_data(
        {
            "id": "KV_id",
            "name": "abc_xyz",
            "key": "keyX",
            "value": "45",
        }
    )
    assert isinstance(keyvaluedict, KeyValueEntry)
    assert keyvaluedict.id == "KV_id"
    assert keyvaluedict.name == "abc_xyz"
    assert keyvaluedict.key == "keyX"
    assert keyvaluedict.value == "45"
    assert keyvaluedict.uri is None


def test_keyvaluepairs():
    entry1 = {
        "id": "id1",
        "name": "abc_xyz",
        "key": "key1",
        "value": "45",
        "uri": "https://tinyurl.com/ybyr7b62",
    }
    entry2 = KeyValueEntry.from_data(
        {
            "id": "id2",
            "name": "abc_xyz2",
            "key": "key2",
            "value": "48",
            "uri": EXAMPLE_URL,
        }
    )

    kvp = KeyValuePairs(entries=[entry1, entry2])
    assert len(kvp) == 2
    for key in ["key1", "key2"]:
        assert key in kvp
    assert kvp["key2"] == entry2
