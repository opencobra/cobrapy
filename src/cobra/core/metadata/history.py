"""Encodes History and Creator.

The history allows to encode provenance meta-data about
model objects. The history allows to encode who created or modified
objects in a model with respective time stamps.
"""
from datetime import datetime
from typing import Dict, Iterable, List, NamedTuple, Optional, Union


STRTIME_FORMAT = "%Y-%m-%dT%H:%M:%S%z"


class History:
    """History object encoding object provenance.

    Parameters
    ----------
    creators : list
        A list of Creators
    created_date : datetime
        The datetime of creation in W3CDTF ISO 8601 format
    modified_dates : list
        A list of datetimes when the object was modified.
    """

    def __init__(
        self,
        creators: Optional[List["Creator"]] = None,
        created_date: Optional[Union[datetime, str]] = None,
        modified_dates: Optional[List[Union[datetime, str]]] = None,
    ):
        """Initialize the class.

        Parameters
        ----------
        creators: list
            list of Creator class. Optional, default None.
        created_date: datetime
            Created date, in HistoryDateTime class. Optional, default None.
        modified_dates: list
            Dates when this annotation was modified. List of datetime dates or strings.
            Optional, default None.
        """
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
    def creators(self) -> List["Creator"]:
        """Get creators for History.

        Returns
        -------
        list
            A list of Creator objects.
        """
        return self._creators

    @creators.setter
    def creators(self, values: Iterable[Union[Dict, "Creator"]]) -> None:
        """Set creators for History.

        Parameters
        ----------
        values: iterable
            An iterable of dictionaries and/or Creator objects.
        """
        self._creators = [Creator.from_data(v) for v in values]

    @staticmethod
    def parse_datetime(value: Optional[Union[str, datetime]]) -> Optional[datetime]:
        """Parse datetime into str.

        Parameters
        ----------
        value: str or datetime
            Optional. If None, the function will return None.
            str is converted to datetime.
            If given datetime, the format will be validated.

        Returns
        -------
        datetime: optional
            Returns None if given None.

        Raises
        ------
        TypeError
            If value is not None, or an instance of str, datetime.
        """
        if value is None:
            return None
        if isinstance(value, datetime):
            return value
        elif isinstance(value, str):
            return History.date_from_str(value)
        else:
            raise TypeError(
                f"Invalid type passed for datetime. "
                f"Accepted types are 'str' or 'datetime' objects: {value}"
            )


    @staticmethod
    def date_from_str(datetime_str: str) -> datetime:
        """Validate if the date format is of type w3cdtf ISO 8601.

        Parameters
        ----------
        datetime_str: str
            Datetime in string format.

        Returns
        -------
        datetime if valid format

        Raises
        ------
        ValueError if not valid.
        """
        if not isinstance(datetime_str, str):
            raise TypeError(
                f"The date passed must be of type string: {datetime_str}")

        # python 3.6 doesn't allow : (colon) in the utc offset.
        try:
            return datetime.strptime(datetime_str, STRTIME_FORMAT)
        except ValueError as e:
            # checking for python 3.6
            if "Z" in datetime_str:
                try:
                    return datetime.strptime(
                        datetime_str.replace("Z", ""), "%Y-%m-%dT%H:%M:%S"
                    )
                except ValueError as e1:
                    raise ValueError(str(e1))
            else:
                utcoff = datetime_str[20:25]
                utcoff_p36 = utcoff.replace(":", "")
                date_p36 = datetime_str.replace(utcoff, utcoff_p36)
                try:
                    return datetime.strptime(date_p36, STRTIME_FORMAT)
                except ValueError:
                    raise ValueError(str(e))

    @property
    def created_date(self) -> Optional[datetime]:
        """Get created date for History.

        Returns
        -------
        datetime
        """
        return self._created_date

    @created_date.setter
    def created_date(self, date: Union[str, "datetime"]) -> None:
        """Set created date for History.

        Parameters
        ----------
        date: str or HistoryDateTime
        """
        self._created_date = self.parse_datetime(date)

    @property
    def modified_dates(self) -> List[datetime]:
        """Get modified dates.

        Returns
        -------
        list
            List of datetimes when this annotation was modified, if any exist.
            List can be empty.
        """
        return self._modified_dates

    @modified_dates.setter
    def modified_dates(self, dates: Iterable[Union[str, datetime]]) -> None:
        """Set modified dates.

        Parameters
        -------
        list
            List of HistoryDateTimes or dictionaries when this annotation was modified.
        """
        self._modified_dates = [self.parse_datetime(d) for d in dates]

    @staticmethod
    def from_data(data: Union[Dict, "History"]) -> "History":
        """Parse history from data.

        Parameters
        ----------
        data: dict or History
            Dict will be parsed to History object.

        Returns
        -------
        History

        Raises
        ------
        TypeError
            If data is neither dict, History or None.
        """
        if data is None:
            return History()
        elif isinstance(data, History):
            return data
        elif isinstance(data, dict):
            return History(**data)
        else:
            raise TypeError(f"Unsupported type for History: '{data}'")

    def is_empty(self) -> bool:
        """Check if history is empty.

        Returns
        -------
        bool
            Returns False if at least one history attribute is set, else True.
        """
        if self.creators:
            return False
        if self.created_date:
            return False
        if self.modified_dates:
            return False
        return True

    def __eq__(self, history: "History") -> bool:
        """Check equality of two history objects.

        A history is equal if all attributes are equal.

        Returns
        -------
        bool - True if equal, False otherwise.
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

    def to_dict(self) -> Dict:
        """Return dictionary representation of History.

        Returns
        -------
        dict - Dictionary representation, of this format
        {
        "creators": list[dict]
        "created_date": str
        "modified_dates": list[str]
        }
        """
        return {
            "creators": [c.to_dict() for c in self.creators],
            "created_date": self.created_date.isoformat(),
            "modified_dates": [mod_date.isoformat() for mod_date in self._modified_dates],
        }

    def __str__(self) -> str:
        """Return a string representation.

        Returns
        -------
        str
            History in a flattened dictionary.
        """
        return str(self.to_dict())

    def __repr__(self):
        """Return a string with module and class name.

        Returns
        -------
        str
            History in a string, with module and class name.
        """
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({self.creators}, {self.created_date}, {self.modified_dates})"
        )


class Creator(NamedTuple):
    """Metadata for person who created an object.

    Parameters
    ----------
    given_name: str
        Optional. Default None.
    family_name: str
        Optional. Default None.
    email: str
        Optional. Default None.
    organisation: str
        Optional. Default None.
    """

    given_name: Optional[str] = None
    family_name: Optional[str] = None
    email: Optional[str] = None
    organisation: Optional[str] = None

    @staticmethod
    def from_data(data: Union[Dict, "Creator"]) -> "Creator":
        """Parse creator from data.

        Parameters
        ----------
        data: dict or Creator
            Dictionary will be converted to Creator class.

        Returns
        -------
        Creator - the creator in the Creator class
        """
        if not data:
            return Creator()
        elif isinstance(data, Creator):
            return data
        elif isinstance(data, dict):
            return Creator._make(
                [
                    data["given_name"],
                    data["family_name"],
                    data["email"],
                    data["organisation"],
                ]
            )
        else:
            raise TypeError(f"Invalid format for Creator: {data}")

    def to_dict(self) -> Dict:
        """Convert Creator to dictionary.

        Returns
        -------
        dict in this format
        {
            "given_name": str,
            "family_name": str,
            "email": str,
            "organisation": str,
        }
        """
        return dict(self._asdict())

    def __str__(self) -> str:
        """Return string representation of Creator.

        Returns
        -------
        str
            String version of flattened dictionary.
        """
        return str(self.to_dict())

    def __repr__(self):
        """Return the Creator with module, class, and internal fields.

        Returns
        -------
        str
        """
        return (
            f"{self.__class__.__module__}.{self.__class__.__qualname__}"
            f"({self.given_name} {self.family_name} {self.email} "
            f"{self.organisation})"
        )
