from datetime import datetime
from typing import Dict, Iterable, List, Union


class HistoryDateTime:
    """ Class representation of DateTimes allowed inside model history.
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
    def datetime(self) -> str:
        return self._datetime

    @datetime.setter
    def datetime(self, value: str) -> None:
        """ Before setting the date, checks if date is in valid format or not. """
        if value is None:
            self._datetime = None
        elif isinstance(value, str):
            self.validate_date(value)
            self._datetime = value
        elif isinstance(value, datetime):
            self._datetime = value.strftime("%Y-%m-%dT%H:%M:%S%z")
        else:
            raise TypeError(
                f"Invalid type passed for datetime. "
                f"Accepted types are 'str' or "
                f"'datetime' objects: {value}"
            )

    def set_current_datetime(self) -> None:
        utcoffset = "+00:00"
        current_datetime = datetime.utcnow()
        self._datetime = current_datetime.strftime("%Y-%m-%dT%H:%M:%S%z")
        self._datetime += utcoffset

    @staticmethod
    def validate_date(date_text: str) -> bool:
        """Validate if the date format is of type w3cdtf ISO 8601"""
        if not isinstance(date_text, str):
            raise TypeError(f"The date passed must be of " f"type string: {date_text}")

        # python 3.6 doesn't allow : (colon) in the utc offset.
        try:
            datetime.strptime(date_text, "%Y-%m-%dT%H:%M:%S%z")
        except ValueError as e:
            # checking for python 3.6
            if 'Z' in date_text:
                try:
                    datetime.strptime(date_text.replace("Z", ""), "%Y-%m-%dT%H:%M:%S")
                except ValueError as e1:
                    raise ValueError(str(e1))
                return True
            else:
                utcoff = date_text[20:25]
                utcoff_p36 = utcoff.replace(":", "")
                date_p36 = date_text.replace(utcoff, utcoff_p36)
                try:
                    datetime.strptime(date_p36, "%Y-%m-%dT%H:%M:%S%z")
                except ValueError as e2:
                    raise ValueError(str(e2))
                return True

        return True

    def __eq__(self, historydatetime_obj: "HistoryDateTime") -> bool:
        return self.datetime == historydatetime_obj.datetime

    def __str__(self) -> str:
        return str(self.datetime)

    def __repr__(self) -> str:
        return self.__str__()


class History:
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

    def __init__(
        self,
        creators: List = None,
        created_date: HistoryDateTime = None,
        modified_dates: List = None,
    ):
        self._modified_dates = []

        self.creators = creators
        self.created_date = created_date
        self.modified_dates = modified_dates

    @staticmethod
    def parse_history(data: Union[Dict, "History"]) -> "History":
        if data is None:
            return History()
        elif isinstance(data, History):
            return data
        elif isinstance(data, dict):
            for key in ["creators", "created_date", "modified_dates"]:
                if key not in data:
                    data[key] = None
            return History(**data)
        else:
            raise TypeError(f"Invalid format for History: '{data}'")

    @property
    def creators(self) -> List:
        return self._creators

    @creators.setter
    def creators(self, value: List) -> None:
        if value is None:
            value = []
        if not isinstance(value, list):
            raise TypeError(f"'creators' must be a list type: {value}")
        else:
            self._creators = value

    @property
    def created_date(self) -> "HistoryDateTime":
        return self._created_date

    @created_date.setter
    def created_date(self, value: Union[str, HistoryDateTime]) -> None:
        if value is None:
            self._created_date = None
        elif isinstance(value, str):
            self._created_date = HistoryDateTime(value)
        elif isinstance(value, HistoryDateTime):
            self._created_date = value
        else:
            raise TypeError(f"'created_date' must be of type DateTime: {value}")

    @property
    def modified_dates(self) -> List:
        return self._modified_dates

    @modified_dates.setter
    def modified_dates(self, value: Iterable) -> None:
        if value is None:
            value = []

        if not isinstance(value, Iterable):
            raise TypeError(f"'modified_dates' must be an Iterable: {value}")

        for date in value:
            if isinstance(date, str):
                self._modified_dates.append(HistoryDateTime(date))
            elif isinstance(date, HistoryDateTime):
                self._modified_dates.append(date)
            else:
                raise TypeError(f"'modified_date' must be of type DateTime: {date}")

    def is_set_history(self) -> bool:
        """Checks if history is set.
        Returns true if at least one attribute is set.
        """
        return len(self.creators) == 0 and \
            self.created_date is not None and \
            len(self.modified_dates) != 0

    def __eq__(self, history_obj: "History") -> bool:
        """ Checking equality of two history objects. """
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

    def __str__(self) -> str:
        return str(
            {
                "creators": self.creators,
                "created_date": self.created_date.datetime,
                "modified_dates": [
                    modified_date.datetime for modified_date in self.modified_dates
                ],
            }
        )

    def __repr__(self) -> str:
        return self.__str__()


class Creator:
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

    def __init__(
        self,
        first_name: str = None,
        last_name: str = None,
        email: str = None,
        organization_name: str = None,
    ):

        self.first_name = first_name
        self.last_name = last_name
        self.email = email
        self.organization_name = organization_name

    @staticmethod
    def parse_creator(data: Union[Dict, "Creator"]) -> "Creator":
        if data is None:
            return Creator()
        elif isinstance(data, dict):
            return Creator(**data)
        elif isinstance(data, Creator):
            return data
        else:
            raise TypeError(f"Invalid format for Creator: {data}")

    def __eq__(self, creator_obj: "Creator") -> bool:
        """ Compare data of two Creator object. If anything
        if found unequal, return False, else return True.
        """
        if (
            not self.first_name == creator_obj.first_name
            or not self.last_name == creator_obj.last_name
            or not self.email == creator_obj.email
            or not self.organization_name == creator_obj.organization_name
        ):
            return False
        else:
            return True

    def __str__(self) -> str:
        return str(
            {
                "first_name": self.first_name,
                "last_name": self.last_name,
                "email": self.email,
                "organization_name": self.organization_name,
            }
        )

    def __repr__(self) -> str:
        return self.__str__()
