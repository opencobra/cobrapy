# -*- coding: utf-8 -*-

from __future__ import absolute_import

import datetime
from typing import List


class History(object):
    """
    Class representation of history of a given component.
    The history allows to store creator, created date and
    modification dates.

    Parameters
    ----------
    creators : list
        A list of cobra creators
    created : string
        The date when component is created in W3CDTF ISO 8601 format
    modified : list
        A list of dates (of type DateTime) about the component modification

    Attributes
    ----------
    creator : list
        A list containing details of creators of the component
    created : DateTime
        The date when component is created in W3CDTF ISO 8601 format
    modified : list
        A list of dates about the component modification

    """

    def __init__(self, creators: 'list' = None, created: 'DateTime' = None,
                 modified: 'list' = None):
        if creators is None:
            creators = []
        self._creators = list(creators)
        if modified is None:
            modified = []
        self._modified = list(modified)
        if created is None:
            self._created = None
        else:
            self._created = created

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
    def created(self):
        return getattr(self, "_created", None)

    @created.setter
    def created(self, value):
        if not isinstance(value, DateTime):
            raise TypeError("Created date must be of type"
                            " DateTime: {}".format(value))
        else:
            self._created = value

    @property
    def creators(self):
        return getattr(self, "_creators", [])

    @creators.setter
    def creators(self, value):
        if not isinstance(value, list):
            raise TypeError("Creators must be wrapped inside "
                            "a list: {}".format(value))
        else:
            self._creators = value

    @property
    def modified(self):
        return getattr(self, "_modified", [])

    @modified.setter
    def modified(self, value):
        if not isinstance(value, list):
            raise TypeError("Modified dates must be wrapped inside"
                            " a list: {}".format(value))
        else:
            self._modified = value

    def isSetHistory(self):
        if self.created is None and len(self.creators) == 0 \
           and len(self.modified) == 0:
                return False
        else:
            return True

    def equals(self, history):
        # checking equality of creators
        if len(self.creators) != len(history.creators):
            return False
        creator_len = len(self.creators)
        for index in range(creator_len):
            if not self.creators[index].equals(history.creators[index]):
                return False

        # checking equality of created date
        if not self.created.equals(history.created):
            return False

        # checking equality of modified
        if len(self.modified) != len(history.modified):
            return False
        modified_len = len(self.modified)
        for index in range(modified_len):
            if not self.modified[index].equals(history.modified[index]):
                return False

        return True

    def __str__(self):
        return str({
            "creators": self.creators,
            "created": self.created.getDateString(),
            "modified": [(modified_date.getDateString()) for modified_date
                         in self.modified]})

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

    def equals(self, creator):
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


class DateTime(object):
    """
    Class representation of dates allowed inside model history.
    This class make sure that dates passed must be of the form :
    %Y-%m-%dT%H:%M:%S%z

    Parameter
    ---------
    date_text : str
        date in the form of a string
    """

    def __init__(self, date_text: 'str' = None):
        if date_text is None:
            date_text = "2000-01-01T00:00:00+00:00"
        self.validateDate(date_text)
        self._date = date_text

    def getDateString(self):
        return getattr(self, "_date", None)

    def setDateFromString(self, value):
        """
        Before setting the date, it first checks if date
        is in valid format or not.
        """
        self.validateDate(value)
        self._date = value

    def equals(self, datetime):
        if self.getDateString() == datetime.getDateString():
            return True
        else:
            return False

    @staticmethod
    def validateDate(date_text):
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
