# -*- coding: utf-8 -*-

import os
from pathlib import Path

import pytest

from cobra.core.metadata.history import Creator, DateTime, History
from cobra.io import read_sbml_model


def _read_ecoli_annotation_model(data_directory):
    test_xml = os.path.join(data_directory, "e_coli_core_for_annotation.xml")
    model = read_sbml_model(test_xml)
    return model


def test_create_history():
    history = History(
        creators=[
            Creator(
                first_name="Matthias",
                last_name="Koenig",
                organization_name="HU",
                email="test@test.com"
            ),
            Creator(
                first_name="Andreas",
                last_name="Draeger",
                organization_name="University of Tübingen",
                email="test2@test2.com"
            )
        ],
        created=DateTime("2020-06-26T02:34:30+05:30"),
        modified=[
            DateTime("2020-06-26T12:34:11+00:00"),
            DateTime("2020-06-26T00:34:11+05:30")
        ]
    )
    assert len(history.creators) == 2
    assert type(history.created) == DateTime
    assert history.created.getDateString() == "2020-06-26T02:34:30+05:30"
    assert len(history.modified) == 2


def test_history_from_ecoli_xml(data_directory):
    model = _read_ecoli_annotation_model(data_directory)
    history = History(
       creators=[
          Creator(
             first_name="Matthias",
             last_name="Koenig",
             email="koenigmx@hu-berlin.de",
             organization_name="Humboldt-University Berlin, "
                               "Institute for Theoretical Biology"
          )
       ],
       created=DateTime("2019-03-06T14:40:55Z"),
       modified=[
          DateTime("2019-03-06T14:40:55Z")
       ]
    )
    assert model.annotation.history.equals(history)


def test_create_creator():
    creator = Creator(
        first_name="Matthias",
        last_name="Koenig",
        organization_name="HU",
        email="test@test.com"
    )
    assert creator.first_name == "Matthias"
    assert creator.last_name == "Koenig"
    assert creator.email == "test@test.com"
    assert creator.organization_name == "HU"


def test_DateTime():
    # valid date
    date = DateTime("2020-06-26T02:34:11+05:30")
    assert date.getDateString() == "2020-06-26T02:34:11+05:30"
    # invalid date (seconds > 59)
    with pytest.raises(ValueError):
        date.setDateFromString("2020-06-26T02:34:70+05:30")
    # valid date
    assert date.setDateFromString("2020-06-26T12:34:11+00:00") is None
