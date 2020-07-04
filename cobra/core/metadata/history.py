# -*- coding: utf-8 -*-

from __future__ import absolute_import

from datetime import datetime
import time
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

    def __init__(self, datetime_obj: [str, datetime] = None):
        self._datetime = None
        self.datetime = datetime_obj

    @property
    def datetime(self):
        return self._datetime

    @datetime.setter
    def datetime(self, value):
        """
        Before setting the date, it first checks if date
        is in valid format or not.
        """
        if value is None:
            self._datetime = None
        elif isinstance(value, str):
            self.validate_date(value)
            self._datetime = value
        elif isinstance(value, datetime):
            self._datetime = value.strftime("%Y-%m-%dT%H:%M:%S%z") \
                             + self.get_utc_offset_str()
        else:
            raise TypeError("Invalid type passed for datetime. "
                            "Accepted types are 'str' or "
                            "'datetime' objects: {}")

    def set_current_datetime(self):
        current_datetime = datetime.now()
        self._datetime = current_datetime.strftime("%Y-%m-%dT%H:%M:%S%z") \
                         + self.get_utc_offset_str()

    def get_utc_offset_str(self):
        """
        Returns a UTC offset string of the current time suitable
        for use in the most widely used timestamps (i.e. ISO 8601,
        RFC 3339). For example: 10 hours ahead, 5 hours behind,
        and time is UTC: +10:00, -05:00, +00:00
        """
        # Calculate the UTC time difference in seconds.
        timestamp = time.time()
        time_now = datetime.fromtimestamp(timestamp)
        time_utc = datetime.utcfromtimestamp(timestamp)
        utc_offset_secs = (time_now - time_utc).total_seconds()
        # Flag variable to hold if the current time is behind UTC.
        is_behind_utc = False
        # If the current time is behind UTC convert the offset
        # seconds to a positive value and set the flag variable.
        if utc_offset_secs < 0:
            is_behind_utc = True
            utc_offset_secs *= -1
        # Build a UTC offset string suitable for use in a timestamp.
        if is_behind_utc:
            pos_neg_prefix = "-"
        else:
            pos_neg_prefix = "+"
        utc_offset = time.gmtime(utc_offset_secs)
        utc_offset_fmt = time.strftime("%H:%M", utc_offset)
        utc_offset_str = pos_neg_prefix + utc_offset_fmt
        return utc_offset_str

    @staticmethod
    def validate_date(date_text):
        """Validate if the date format is of type w3cdtf ISO 8601"""
        if not isinstance(date_text, str):
            raise TypeError("The date passed must be of "
                            "type string: {}".format(date_text))

        try:
            datetime.strptime(date_text, '%Y-%m-%dT%H:%M:%S%z')
        except ValueError as e:
            raise ValueError(str(e))
        return True

    def __eq__(self, historydatetime_obj):
        return self.datetime == historydatetime_obj.datetime

    def __str__(self):
        return str(self.datetime)

    def __repr__(self):
        return self.__str__()


class History(object):
    """
    Class representation of history of a given component.
    The history allows to store creator, created date and
    modification dates.

    Parameters
    ----------
    creators : list
        A list of cobra creators
    created_date : HistoryDateTime
        The date when component is created in W3CDTF ISO 8601 format
    modified_dates : list
        A list of dates (of type DateTime) about the component modification

    Attributes
    ----------
    creators : list
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
            for key in ["creators", "created_date",
                        "modified_dates"]:
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
            raise TypeError("'creators' must be a list "
                            "{}".format(value))
        else:
            self._creators = value

    @property
    def created_date(self):
        return self._created_date

    @created_date.setter
    def created_date(self, value):
        if value is None:
            self._created_date = None
        elif isinstance(value, str):
            self._created_date = HistoryDateTime(value)
        elif not isinstance(value, HistoryDateTime):
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
            raise TypeError("'modified_dates' must be "
                            " a list: {}".format(value))
        else:
            for date in value:
                if isinstance(date, str):
                    self._modified_dates.append(HistoryDateTime(date))
                elif isinstance(date, HistoryDateTime):
                    self._modified_dates.append(date)
                else:
                    raise TypeError("'modified_date' must be of type"
                                    " DateTime: {}".format(date))

    def is_set_history(self) -> bool:
        """Checks if history is set.
        Returns true if at least one attribute is set.
        """
        if self.creators or self.created_date or \
                self.modified_dates:
            return True
        return False

    def __eq__(self, history_obj):
        """ Checking equality of two history objects.
        """
        # check creators
        if len(self.creators) != len(history_obj.creators):
            return False
        for k, creator in enumerate(self.creators):
            if not creator == history_obj.creators[k]:
                return False

        # checking equality of created date
        if not self.created_date == history_obj.created_date:
            return False

        # checking equality of modified
        if len(self.modified_dates) != len(history_obj.modified_dates):
            return False
        for k, modified_date in enumerate(self.modified_dates):
            if not modified_date == history_obj.modified_dates[k]:
                return False

        return True

    def __str__(self):
        return str({
            "creators": self.creators,
            "created_date": self.created_date.datetime,
            "modified_dates": [modified_date.datetime for modified_date
                               in self.modified_dates]})

    def __repr__(self):
        return self.__str__()


class Creator(object):
    """Class representation of a Creator

    Parameters
    ----------
        first_name : str,
        last_name : str,
        email : str,
        organization_name : str
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

    def __eq__(self, creator_obj):
        """ Compare data of two Creator object. If anything
        if found unequal, return False, else return True.
        """
        if not self.first_name == creator_obj.first_name or \
           not self.last_name == creator_obj.last_name or \
           not self.email == creator_obj.email or \
           not self.organization_name == creator_obj.organization_name:
            return False
        else:
            return True

    def __str__(self):
        return str({
            "first_name": self.first_name,
            "last_name": self.last_name,
            "email": self.email,
            "organization_name": self.organization_name}
        )

    def __repr__(self):
        return self.__str__()

