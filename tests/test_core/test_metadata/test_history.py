import os
from datetime import datetime

import pytest

from cobra.core.metadata.history import Creator, History, HistoryDatetime
from cobra.io import read_sbml_model


def _read_ecoli_annotation_model(data_directory):
    """Helper function to read model with history elements."""
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
                email="test@test.com",
            ),
            Creator(
                given_name="Andreas",
                family_name="Draeger",
                organisation="University of Tübingen",
                email="test2@test2.com",
            ),
        ],
        created_date=HistoryDatetime("2020-06-26T02:34:30+05:30"),
        modified_dates=[
            HistoryDatetime("2020-06-26T12:34:11+00:00"),
            HistoryDatetime("2020-06-26T00:34:11+05:30"),
        ],
    )
    assert len(history.creators) == 2
    assert isinstance(history.created_date, HistoryDatetime)
    assert history.created_date.datetime == "2020-06-26T02:34:30+05:30"
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
            )
        ],
        created_date=HistoryDatetime("2019-03-06T14:40:55Z"),
        modified_dates=[HistoryDatetime("2019-03-06T14:40:55Z")],
    )
    assert model.annotation.history == history


def test_create_creator():
    creator = Creator(
        given_name="Matthias",
        family_name="König",
        organisation="HU",
        email="test@test.com",
    )
    assert creator is not None
    assert creator.given_name == "Matthias"
    assert creator.family_name == "König"
    assert creator.organisation == "HU"
    assert creator.email == "test@test.com"


def test_historydatetime():
    # valid date
    dt_str1 = "2020-06-26T02:34:11+05:30"
    datetime_obj = HistoryDatetime(dt_str1)
    assert datetime_obj.datetime == dt_str1

    # valid date
    datetime_obj.datetime = "2020-06-26T12:34:11+00:00"
    assert datetime_obj.datetime == "2020-06-26T12:34:11+00:00"
    datetime_obj.datetime = None
    assert datetime_obj.datetime is None

    # create from python datetime
    datetime_obj.datetime = datetime.now()
    assert datetime_obj.datetime is not None
