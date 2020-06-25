from cobra.core.model import Object
from cobra.core.metadata.history import History, Creator


def test_create_history():
    h = History(

    )


    assert 0 == 1

def test_create_creator():
    h = Creator(
        first_name="Matthias",
        last_name="Koenig",
        organization_name="HU",
        email="test@test.com"
    )
    assert h.first_name == "Matthias"
    # FIXME: simple test;