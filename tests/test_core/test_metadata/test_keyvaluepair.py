from cobra.core.metadata.custompairs import KeyValueEntry, KeyValuePairs


def test_keyvalueentry():
    keyvaluedict = KeyValueEntry.from_data(
        {
            "id": "KV_id",
            "name": "abc_xyz",
            "key": "keyX",
            "value": "45",
            "uri": "https://tinyurl.com/ybyr7b62",
        }
    )
    assert isinstance(keyvaluedict, KeyValueEntry)
    assert keyvaluedict.id == "KV_id"
    assert keyvaluedict.name == "abc_xyz"
    assert keyvaluedict.key == "keyX"
    assert keyvaluedict.value == "45"
    assert keyvaluedict.uri == "https://tinyurl.com/ybyr7b62"


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
            "uri": "https://tinyurl2.com/ybyr7b62",
        }
    )

    kvp = KeyValuePairs(entries=[entry1, entry2])
    print(kvp)
    assert len(kvp) == 2
    for key in ["key1", "key2"]:
        print("***", key, "***")
        assert key in kvp
    assert kvp["key2"] == entry2
