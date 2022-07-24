import os
from datetime import datetime, timedelta, timezone

import pytest

from cobra.core.metadata.history import Creator, History
from cobra.io import read_sbml_model


JUNE_26TH = "2020-06-26T12:34:11+00:00"
JUNE_26TH_530 = "2020-06-26T00:34:11+05:30"
TEST_COM = "test@test.com"


def _read_ecoli_annotation_model(data_directory):
    """Read model with history elements."""
    test_xml = os.path.join(data_directory, "e_coli_core_for_annotation.xml")
    model = read_sbml_model(test_xml)
    return model


def test_create_history():
    history = History(
        creators=[
            Creator(
                given_name="Matthias",
                family_name="Koenig",
                organisation="HU",
                email=TEST_COM,
            ),
            Creator(
                given_name="Andreas",
                family_name="Draeger",
                organisation="University of Tübingen",
                email="test2@test2.com",
            ),
        ],
        created_date=JUNE_26TH_530,
        modified_dates=[
            JUNE_26TH,
            JUNE_26TH_530,
        ],
    )
    assert len(history.creators) == 2
    assert history.created_date.isoformat() == JUNE_26TH_530
    assert history.created_date == History.parse_datetime(JUNE_26TH_530)
    assert len(history.modified_dates) == 2


def test_history_from_ecoli_xml(data_directory):
    model = _read_ecoli_annotation_model(data_directory)
    history = History(
        creators=[
            Creator(
                given_name="Matthias",
                family_name="Koenig",
                email="koenigmx@hu-berlin.de",
                organisation="Humboldt-University Berlin, "
                "Institute for Theoretical Biology",
            ),
            Creator(
                given_name="Matthias",
                family_name="Koenig",
                email="koenigmx@hu-berlin.de",
            ),
            Creator(
                given_name="Matthias",
                family_name="Koenig",
                organisation="Humboldt-University Berlin, "
                "Institute for Theoretical Biology",
            ),
            Creator(
                given_name="Matthias",
                family_name="Koenig",
            ),
        ],
        created_date="2019-03-06T14:40:55Z",
        modified_dates=[
            "2019-03-06T14:40:55Z",
            "2019-03-06T14:41:55Z",
        ],
    )
    assert model.annotation.history == history
    model.annotation.history.created_date = None
    assert model.annotation.history == History(
        creators=[
            Creator(
                given_name="Matthias",
                family_name="Koenig",
                email="koenigmx@hu-berlin.de",
                organisation="Humboldt-University Berlin, "
                "Institute for Theoretical Biology",
            ),
            Creator(
                given_name="Matthias",
                family_name="Koenig",
                email="koenigmx@hu-berlin.de",
            ),
            Creator(
                given_name="Matthias",
                family_name="Koenig",
                organisation="Humboldt-University Berlin, "
                "Institute for Theoretical Biology",
            ),
            Creator(
                given_name="Matthias",
                family_name="Koenig",
            ),
        ],
        modified_dates=[
            "2019-03-06T14:40:55Z",
            "2019-03-06T14:41:55Z",
        ],
    )


def test_create_creator():
    creator = Creator(
        given_name="Matthias",
        family_name="König",
        organisation="HU",
        email=TEST_COM,
    )
    assert creator.given_name == "Matthias"
    assert creator.family_name == "König"
    assert creator.organisation == "HU"
    assert creator.email == TEST_COM

    creator = Creator(
        given_name="Matthias",
        organisation="HU",
        email=TEST_COM,
    )
    assert creator.given_name == "Matthias"
    assert creator.family_name is None
    assert creator.organisation == "HU"
    assert creator.email == TEST_COM

    creator = Creator(
        **{
            "given_name": "Matthias",
            "family_name": "König",
            "organisation": "HU",
            "email": TEST_COM,
        }
    )

    assert creator.given_name == "Matthias"
    assert creator.family_name == "König"
    assert creator.organisation == "HU"
    assert creator.email == TEST_COM

    creator = Creator(
        **{
            "given_name": "Matthias",
            "organisation": "HU",
            "email": TEST_COM,
        }
    )

    assert creator.given_name == "Matthias"
    assert creator.family_name is None
    assert creator.organisation == "HU"
    assert creator.email == TEST_COM

    creator = Creator().from_data(
        {
            "given_name": "Matthias",
            "family_name": "König",
            "organisation": "HU",
            "email": TEST_COM,
        }
    )

    assert creator.given_name == "Matthias"
    assert creator.family_name == "König"
    assert creator.organisation == "HU"
    assert creator.email == TEST_COM

    creator = Creator().from_data(
        {
            "given_name": "Matthias",
            "organisation": "HU",
            "email": TEST_COM,
        }
    )

    assert creator.given_name == "Matthias"
    assert creator.family_name is None
    assert creator.organisation == "HU"
    assert creator.email == TEST_COM


def test_historydatetime():
    # valid date
    dt_str1 = JUNE_26TH_530
    datetime_obj = History(created_date=dt_str1)
    assert datetime_obj.created_date.isoformat() == dt_str1

    # valid date
    datetime_obj = History(created_date=JUNE_26TH)
    assert datetime_obj.created_date.isoformat() == JUNE_26TH
    datetime_obj.created_date = None
    assert datetime_obj.created_date is None

    # create from python datetime
    datetime_obj.created_date = datetime.now()
    assert datetime_obj.created_date is not None

    # incorrect format
    with pytest.raises(ValueError):
        datetime_obj.created_date = datetime.now(timezone(timedelta(hours=1))).strftime(
            "%Y-%m-%dT%H:%M:%S%Z"
        )

    with pytest.raises(TypeError):
        datetime_obj.created_date = {"date": "June 26th"}
