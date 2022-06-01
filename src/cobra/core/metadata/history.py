"""Encodes History and Creator.

The history allows to encode provenance meta-data about
model objects. The history allows to encode who created or modified
objects in a model with respective time stamps.
"""
from datetime import datetime
from typing import Dict, Iterable, List, Optional, Union


STRTIME_FORMAT = "%Y-%m-%dT%H:%M:%S%z"


class History:
    """History object encoding object provenance.

    Parameters
    ----------
    creators : list
        A list of Creators
    created_date : HistoryDatetime
        The datetime of creation in W3CDTF ISO 8601 format
    modified_dates : list
        A list of datetimes when the object was modified.
    """

    def __init__(
        self,
        creators: List["Creator"] = None,
        created_date: "HistoryDatetime" = None,
        modified_dates: List["HistoryDatetime"] = None,
    ):
        if modified_dates is None:
            modified_dates = []
        if creators is None:
            creators = []
        self._creators = []
        self._created_date = None
        self._modified_dates = []

        # use properties to set fields
        self.creators = creators
        self.created_date = created_date
        self.modified_dates = modified_dates

    @property
    def creators(self) -> List:
        return self._creators

    @creators.setter
    def creators(self, values: Iterable["Creator"]) -> None:
        self._creators = [Creator.from_data(v) for v in values]

    @property
    def created_date(self) -> "HistoryDatetime":
        return self._created_date

    @created_date.setter
    def created_date(self, date: Union[str, "HistoryDateTime"]) -> None:
        self._created_date = HistoryDatetime(date)

    @property
    def modified_dates(self) -> List:
        return self._modified_dates

    @modified_dates.setter
    def modified_dates(self, dates: Iterable[Union[str, "HistoryDateTime"]]) -> None:
        self._modified_dates = [HistoryDatetime(d) for d in dates]

    @staticmethod
    def from_data(data: Union[Dict, "History"]) -> "History":
        """Parse history from data."""
        if data is None:
            return History()
        elif isinstance(data, History):
            return data
        elif isinstance(data, dict):
            return History(**data)
        else:
            raise TypeError(f"Unsupported type for History: '{data}'")

    def is_empty(self) -> bool:
        """Checks if history is empty.

        Returns False if at least one history attribute is set, else True.
        """
        if self.creators:
            return False
        if self.created_date.datetime:
            return False
        if self.modified_dates:
            return False
        return True

    def __eq__(self, history: "History") -> bool:
        """Checking equality of two history objects.

        A history is equal if all attributes are equal.
        """
        # check equality of creators
        if len(self.creators) != len(history.creators):
            return False
        for k, creator in enumerate(self.creators):
            if creator != history.creators[k]:
                return False

        # checking equality of created_date
        if self.created_date != history.created_date:
            return False

        # checking equality of modified_dates
        if len(self.modified_dates) != len(history.modified_dates):
            return False
        for k, modified_date in enumerate(self.modified_dates):
            if modified_date != history.modified_dates[k]:
                return False

        return True

    def to_dict(self):
        """Returns dictionary representation."""
        modified_dates = []
        for modified_date in self._modified_dates:
            modified_dates.append(modified_date.datetime)
        return {
            "creators": [c.to_dict() for c in self.creators],
            "created_date": self.created_date.datetime,
            "modified_dates": modified_dates,
        }

    def __str__(self) -> str:
        return str(self.to_dict())

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.creators}>"


class Creator:
    """Metadata for person who created an object.

    Parameters
    ----------
        given_name : str,
        family_name : str,
        email : str,
        organisation : str
    """

    def __init__(
        self,
        given_name: str = None,
        family_name: str = None,
        email: str = None,
        organisation: str = None,
    ):
        self.given_name = given_name  # type: str
        self.family_name = family_name  # type: str
        self.email = email  # type: str
        self.organisation = organisation  # type: str

    @staticmethod
    def from_data(data: Union[Dict, "Creator"]) -> "Creator":
        """Parse creator from data."""
        if not data:
            return Creator()
        elif isinstance(data, Creator):
            return data
        elif isinstance(data, dict):
            return Creator(**data)
        else:
            raise TypeError(f"Invalid format for Creator: {data}")

    def __eq__(self, creator_obj: "Creator") -> bool:
        """Compare creator objects for equality

        Two creators are equal if all fields are equal.
        """
        if (
            self.given_name != creator_obj.given_name
            or self.family_name != creator_obj.family_name
            or self.email != creator_obj.email
            or self.organisation != creator_obj.organisation
        ):
            return False

        return True

    def to_dict(self):
        return {
            "given_name": self.given_name,
            "family_name": self.family_name,
            "email": self.email,
            "organisation": self.organisation,
        }

    def __str__(self) -> str:
        return str(self.to_dict())

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.given_name} {self.family_name}>"


class HistoryDatetime:
    """Datetime allowed in a model history.

    This class make sure that datetimes are of the form:
    %Y-%m-%dT%H:%M:%S%z

    Parameter
    ---------
    datetime: str, datetime
        date in the form of a string
    """

    def __init__(self, history_datetime: str = None):
        self._datetime: Optional[str] = None
        self.datetime = history_datetime

    @property
    def datetime(self) -> str:
        return self._datetime

    @datetime.setter
    def datetime(self, value: str) -> None:
        self._datetime = self.parse_datetime(value)

    def parse_datetime(self, value: str) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, HistoryDatetime):
            return value.datetime
        elif isinstance(value, str):
            self.validate_datetime(value)
            return value
        elif isinstance(value, datetime):
            return value.strftime(STRTIME_FORMAT)
        else:
            raise TypeError(
                f"Invalid type passed for datetime. "
                f"Accepted types are 'str' or "
                f"'datetime' objects: {value}"
            )

    @staticmethod
    def utcnow() -> "HistoryDatetime":
        """HistoryDatetime with current UTC time."""
        return HistoryDatetime(datetime.utcnow().strftime(STRTIME_FORMAT))

    @staticmethod
    def validate_datetime(datetime_str: str) -> None:
        """Validate if the date format is of type w3cdtf ISO 8601.

        Raises ValueError if not valid.
        """
        if not isinstance(datetime_str, str):
            raise TypeError(f"The date passed must be of type string: {datetime_str}")

        # python 3.6 doesn't allow : (colon) in the utc offset.
        try:
            datetime.strptime(datetime_str, STRTIME_FORMAT)
        except ValueError as e:
            # checking for python 3.6
            if "Z" in datetime_str:
                try:
                    datetime.strptime(
                        datetime_str.replace("Z", ""), "%Y-%m-%dT%H:%M:%S"
                    )
                except ValueError as e1:
                    raise ValueError(str(e1))
            else:
                utcoff = datetime_str[20:25]
                utcoff_p36 = utcoff.replace(":", "")
                date_p36 = datetime_str.replace(utcoff, utcoff_p36)
                try:
                    datetime.strptime(date_p36, STRTIME_FORMAT)
                except ValueError:
                    raise ValueError(str(e))

    def __eq__(self, history_datetime: "HistoryDatetime") -> bool:
        return self.datetime == history_datetime.datetime

    def __str__(self) -> str:
        return self.datetime

    def __repr__(self):
        return f"<{self.__class__.__name__} {self.datetime}>"