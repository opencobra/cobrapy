
from cobra.core.species import Species

def test_annotation():
    s = Species()
    print(s.annotation)
    s.annotation["chebi"] = ["1234", "23423432"]
    s.annotation["sbo"] = ["SBO123"]
    print(s.annotation)

    assert "chebi" in s.annotation
    assert "sbo" in s.annotation
    assert len(s.annotation) == 2
    for key in ["keys", "items", "values"]:
        assert hasattr(s.annotation, key)

    # assert 0 == 1

