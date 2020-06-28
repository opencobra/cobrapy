# -*- coding: utf-8 -*-

import datetime
from typing import List


class History(object):
    """
    Class representation of history of a given component.
    The history allows to store creator, created date and
    modification dates.

    Parameters
    ----------
    creator : dict
        A dictionary containing details of creator's name, email and
        organisation name
    created : string
        The date when component is created in W3CDTF ISO 8601 format
    modified : list
        A list of dates about the component modification

    Attributes
    ----------
    creator :Creators
        A dictionary containong details of creator's name, email and
        organisation name
    created : string
        The date when component is created in W3CDTF ISO 8601 format
    modified : list
        A list of dates about the component modification

    """

    def __init__(self, creators: 'list' = None, created: 'DateTime' = None,
                 modified: 'list' = None):
        if creators is None:
            creators = []
        self._creators: List(Creator) = list(creators)
        if modified is None:
            modified = []
        self._modified: List(DateTime) = list(modified)
        if created is None:
            self._created = None
        else:
            self._created = DateTime(created)

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
            raise TypeError("Invalid format for History: '{}'".format(data))

    @property
    def created(self):
        return getattr(self, "_created", None)

    @created.setter
    def created(self, value):
        if not isinstance(value, DateTime):
            raise TypeError("Created date must be of type DateTime: {}".format(value))
        else:
            self._created = value

    @property
    def creators(self):
        return getattr(self, "_creators", [])

    @creators.setter
    def creators(self, value):
        if not isinstance(value, list):
            raise TypeError("Creators must be wrapped inside a list: {}".format(value))
        else:
            self._creators = value

    @property
    def modified(self):
        return getattr(self, "_modified", [])

    @modified.setter
    def modified(self, value):
        if not isinstance(value, list):
            raise TypeError("Modified dates must be wrapped inside a list: {}".format(value))
        else:
            self._modified = value

    def isSetHistory(self):
        if self.created == None and len(self.creators) == 0 and len(self.modified) == 0:
            return False
        else:
            return True

    def __str__(self):
        return str({"creators": self.creators, "created": self.created.getDateString(),
                        "modified": [(modified_date.getDateString()) for modified_date in self.modified]})

    def __repr__(self):
        return self.__str__()

'''
FIXME: remove, use simple list
class ListOfCreators(MutableSequence):
    """A list extension to store each creator's info

    Parameters
    ----------
    creators : list containing info about creators
    """

    def __init__(self, creators=None):
        if creators is None:
            creators = []
        self._sequence = list()
        if not isinstance(creators, list):
            raise TypeError("The data passed for creators must be "
                            "inside a list")
        else:
            for item in creators:
                if isinstance(item, Creator):
                    self._sequence.append(item)
                elif isinstance(item, dict):
                    self._sequence.append(Creator(item))
                else:
                    raise TypeError("Invalid format for Creator: {}".format(item))

    def __len__(self):
        return len(self._sequence)

    def __delitem__(self, index):
        del self._sequence[index]

    def insert(self, index, value):
        if isinstance(value, Creator):
            self._sequence.insert(index, value)
        elif isinstance(value, dict):
            self._sequence.insert(index, Creator(value))
        else:
            raise TypeError("The data passed has invalid format: {}".format(value))

    def append(self, value):
        if isinstance(value, Creator):
            self._sequence.append(value)
        elif isinstance(value, dict):
            self._sequence.append(Creator(value))
        else:
            raise TypeError("The data passed has invalid format: {}".format(value))

    def __setitem__(self, index, value):
        if isinstance(value, Creator):
            self._sequence[index] = value
        elif isinstance(value, dict):
            self._sequence[index] = Creator(value)
        else:
            raise TypeError("The data passed has invalid format: {}".format(value))

    def __getitem__(self, index):
        return self._sequence[index]

    def __str__(self):
        return str(self._sequence)

    def __repr__(self):
        return '{}'.format(self._sequence)
'''

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

    def __init__(self, date_text: 'str' = None):
        if date_text is None:
            date_text = "2000-01-01T00:00:00+00:00"
        self.validateDate(date_text)
        self._date = date_text

    def getDateString(self):
        return self._date

    def setDateFromString(self, value):
        self.validateDate(value)
        self._date = value

    @staticmethod
    def validateDate(date_text):
        """Validate if the date format is of type w3cdtf ISO 8601"""
        if not isinstance(date_text, str):
            raise TypeError("The date passed must be of type string: {}".format(date_text))

        try:
            datetime.datetime.strptime(date_text, '%Y-%m-%dT%H:%M:%S%z')
        except ValueError as e:
            raise ValueError(str(e))
        return True

    def __str__(self):
        return self._date

    def __repr__(self):
        return self._date


# FIXME: remove
# class ModifiedHistory(MutableSequence):
#     """A list extension to store modification dates. Only Restricted
#     type of entries are possible.
#
#     Parameters
#     ----------
#     modifiedList : list containing modification dates in W3CDTF ISO
#                    8601 format
#     """
#
#     def __init__(self, modifiedList=None):
#         if modifiedList is None:
#             modifiedList = []
#         self._sequence = list()
#         if not isinstance(modifiedList, list):
#             raise TypeError("The dates passed must be inside a list: {}".format(modifiedList))
#         for item in modifiedList:
#             if not isinstance(item, str):
#                 raise ValueError("Each date must be of type string: {}".format(item))
#             validateDate(item)
#             self._sequence.append(item)
#
#     def __len__(self):
#         return len(self._sequence)
#
#     def __delitem__(self, index):
#         del self._sequence[index]
#
#     def insert(self, index, value):
#         validateDate(value)
#         self._sequence.insert(index, value)
#
#     def append(self, value):
#         validateDate(value)
#         self._sequence.append(value)
#
#     def __setitem__(self, index, value):
#         validateDate(value)
#         self._sequence[index] = value
#
#     def __getitem__(self, index):
#         return self._sequence[index]
#
#     def __str__(self):
#         return str(self._sequence)
#
#     def __repr__(self):
#         return '{}'.format(self._sequence)
