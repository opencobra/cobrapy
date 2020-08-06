# -*- coding: utf-8 -*-

import pytest

from cobra.core.metadata.keyvaluepair import KeyValueDict, ListOfKeyValue


def test_keyvaluedict():
    keyvaluedict = KeyValueDict.parse_keyValueDict(
        {
            "id": "KV_id",
            "name": "abc_xyz",
            "key": "keyX",
            "value": "45",
            "uri": "https://tinyurl.com/ybyr7b62"
        }
    )
    assert keyvaluedict.id == "KV_id"
    assert keyvaluedict.name == "abc_xyz"
    assert keyvaluedict.key == "keyX"
    assert keyvaluedict.value == "45"
    assert keyvaluedict.uri == "https://tinyurl.com/ybyr7b62"
    # only string type allowed for value
    with pytest.raises(TypeError):
        keyvaluedict.value = 45


def test_listofKeyValue():
    listofkeyvalue = ListOfKeyValue(
        [
            {
                "id": "KV_id",
                "name": "abc_xyz",
                "key": "keyX",
                "value": "45",
                "uri": "https://tinyurl.com/ybyr7b62"
            },
            {
                "id": "KV_id2",
                "name": "abc_xyz2",
                "key": "keyY",
                "value": "48",
                "uri": "https://tinyurl2.com/ybyr7b62"
            }
        ]
    )
    assert len(listofkeyvalue) == 2
