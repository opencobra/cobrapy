# -*- coding: utf-8 -*-

from __future__ import absolute_import

import datetime
from typing import List


class HistoryDateTime(object):
    """
    Class representation of DateTimes allowed inside model history.
    This class make sure that dates passed must be of the form:

    %Y-%m-%dT%H:%M:%S%z

    Parameter
    ---------
    datetime: str
        date in the form of a string
    """

    def __init__(self, datetime=None):
        # FIXME: accept python DateTimes & strings
        self.datetime = datetime

    @property
    def datetime(self):
        return getattr(self, "_date", None)

    @datetime.setter
    def datetime(self, value):
        """
        Before setting the date, it first checks if date
        is in valid format or not.
        """
        # FIXME: add python Datetime support

        if value is None:
            # FIXME: is default a good idea? probably not (better None)
            value = "2000-01-01T00:00:00+00:00"

        self.validate_date(value)
        self._date = value

    def __eq__(self, value):
        return self.datetime == value.datetime

    @staticmethod
    def validate_date(date_text):
        """Validate if the date format is of type w3cdtf ISO 8601"""
        if not isinstance(date_text, str):
            raise TypeError("The date passed must be of "
                            "type string: {}".format(date_text))

        try:
            datetime.datetime.strptime(date_text, '%Y-%m-%dT%H:%M:%S%z')
        except ValueError as e:
            raise ValueError(str(e))
        return True

    def __str__(self):
        return self._date

    def __repr__(self):
        return self._date


class History(object):
    """
    Class representation of history of a given component.
    The history allows to store creator, created date and
    modification dates.

    Parameters
    ----------
    creators : list
        A list of cobra creators
    created_date : string
        The date when component is created in W3CDTF ISO 8601 format
    modified_dates : list
        A list of dates (of type DateTime) about the component modification

    Attributes
    ----------
    creator : list
        A list containing details of creators of the component
    created_date : HistoryDateTime
        The date when component is created in W3CDTF ISO 8601 format
    modified_dates : list
        A list of dates about the component modification

    """

    def __init__(self, creators: List = None, created_date: HistoryDateTime = None,
                 modified_dates: List = None):
        self._creators = []
        self._created_date = None
        self._modified_dates = []

        self.creators = creators
        self.created_date = created_date
        self.modified_dates = modified_dates

    @staticmethod
    def parse_history(data) -> 'History':
        if data is None:
            return History()
        elif isinstance(data, History):
            return data
        elif isinstance(data, dict):
            for key in ["creators", "created", "modified"]:
                if key not in data:
                    data[key] = None
            return History(**data)
        else:
            raise TypeError("Invalid format for "
                            "History: '{}'".format(data))

    @property
    def creators(self):
        return self._creators

    @creators.setter
    def creators(self, value):
        if value is None:
            value = []

        if not isinstance(value, list):
            raise TypeError("Creators must be a list "
                            "{}".format(value))
        else:
            self._creators = value

    @property
    def created_date(self):
        return self._created_date

    @created_date.setter
    def created_date(self, value):
        if not isinstance(value, HistoryDateTime):
            raise TypeError("'created_date' must be of type"
                            " DateTime: {}".format(value))
        else:
            self._created_date = value

    @property
    def modified_dates(self):
        return self._modified_dates

    @modified_dates.setter
    def modified_dates(self, value):
        if value is None:
            value = []
        if not isinstance(value, list):
            raise TypeError("'modified_dates' must be wrapped inside"
                            " a list: {}".format(value))
        else:
            self._modified_dates = value

    def is_set_history(self) -> bool:
        """Checks if history is set.
        Returns true if at least one attribute is set.
        :return:
        """
        if self.creators or self.created_date or self.modified_dates:
            return True
        return False

    def __eq__(self, history):
        """ Checking equality of history elements.

        :param history:
        :return:
        """
        # check creators
        if len(self.creators) != len(history.creators):
            return False
        for k, creator in enumerate(self.creators):
            if not creator == history.creators[k]:
                return False

        # checking equality of created date
        if not self.created_date == history.created_date:
            return False

        # checking equality of modified
        if len(self.modified) != len(history.modified):
            return False
        for k, modified_date in enumerate(self.modified_dates):
            if not modified_date == history.modified_date[k]:
                return False

        return True

    def __str__(self):
        return str({
            "creators": self.creators,
            "created_date": self.created_date.getDateString(),
            "modified_dates": [(modified_date.getDateString()) for modified_date
                         in self.modified_dates]})

    def __repr__(self):
        return self.__str__()


class Creator(object):
    """Class representation of a Creator

    Parameters
    ----------

    Attributes
    ----------
        first_name : str,
        last_name : str,
        email : str,
        organization_name : str
    """
    def __init__(self, first_name: str = None, last_name: str = None,
                 email: str = None, organization_name: str = None):

        self.first_name = first_name
        self.last_name = last_name
        self.email = email
        self.organization_name = organization_name

    @staticmethod
    def parse_creator(data) -> 'Creator':
        if data is None:
            return Creator()
        elif isinstance(data, dict):
            return Creator(**data)
        elif isinstance(data, Creator):
            return data
        else:
            raise TypeError("Invalid format for Creator: {}".format(data))

    def __eq__(self, creator):
        # FIXME: better test via negative or
        if self.first_name == creator.first_name and \
           self.last_name == creator.last_name and \
           self.email == creator.email and \
           self.organization_name == creator.organization_name:
                return True
        else:
            return False

    def __str__(self):
        return str({
            "first_name": self.first_name,
            "last_name": self.last_name,
            "email": self.email,
            "organization_name": self.organization_name}
        )

    def __repr__(self):
        return self.__str__()

