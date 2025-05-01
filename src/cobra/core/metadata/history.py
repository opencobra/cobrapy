"""Encodes History and Creator.

The history allows to encode provenance meta-data about
model objects. The history allows to encode who created or modified
objects in a model with respective time stamps.
"""

import re
from datetime import datetime
from typing import Any, Dict, Iterable, List, Optional, Union


STRTIME_FORMAT = "%Y-%m-%dT%H:%M:%S%z"


class History:
    """History object encoding object provenance.

    Parameters
    ----------
    creators: list
        list of Creator class. Optional, default None.
    created_date: datetime
        Created date. Optional, default None.
    modified_dates: list
        Dates when this annotation was modified. List of datetime dates or strings.
        Optional, default None.
    """

    def __init__(
        self,
        creators: Optional[List["Creator"]] = None,
        created_date: Optional[Union[datetime, str]] = None,
        modified_dates: Optional[List[Union[datetime, str]]] = None,
    ):
        """Initialize the class."""

        self._creators: List[Creator] = []
        self._created_date: Optional[datetime] = None
        self._modified_dates: List[datetime] = []

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
    def creators(self, values: Optional[Iterable[Union[Dict, "Creator"]]]) -> None:
        """Set creators for History.

        Parameters
        ----------
        values: iterable
            An iterable of dictionaries and/or Creator objects.
        """
        if values is None:
            self._creators = []
        else:
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
            raise TypeError(f"The date passed must be of type string: {datetime_str}")

        # python 3.6 doesn't allow : (colon) in the utc offset.
        try:
            datetime_return = datetime.strptime(datetime_str, STRTIME_FORMAT)
            return datetime_return
        except ValueError as e:
            datetime_str = datetime_str.replace("Z", "+0000")
            datetime_str = re.sub(r"(\+\d\d):(\d\d)\Z", "\\1\\2", datetime_str)
            try:
                datetime_return = datetime.strptime(datetime_str, STRTIME_FORMAT)
            except ValueError:
                raise ValueError(str(e))
        return datetime_return

    @property
    def created_date(self) -> Optional[datetime]:
        """Get created date for History.

        Returns
        -------
        datetime
        """
        return self._created_date

    @created_date.setter
    def created_date(self, date: Optional[Union[str, "datetime"]]) -> None:
        """Set created date for History.

        Parameters
        ----------
        date: str or datetime
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
    def modified_dates(self, dates: Optional[Iterable[Union[str, datetime]]]) -> None:
        """Set modified dates.

        Parameters
        -------
        list
            List of datetimes or strings when this annotation was modified.
        """
        if dates is None:
            self._modified_dates = []
        else:
            mds = [self.parse_datetime(d) for d in dates]
            self._modified_dates = [md for md in mds if md is not None]

    @staticmethod
    def from_data(data: Optional[Union[Dict, "History"]]) -> "History":
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

    def __eq__(self, other: "History") -> bool:
        """Check equality of two history objects.

        A history is equal if all attributes are equal.
        If one and only one of self or other is empty will return False.
        If both are empty, will return True.

        Returns
        -------
        bool - True if equal, False otherwise.
        """
        if self.is_empty() and other.is_empty():
            return True
        elif (self.is_empty() and not other.is_empty()) or (
            not self.is_empty() and other.is_empty()
        ):
            return False
        # check equality of creators
        if len(self.creators) != len(other.creators):
            return False
        for k, creator in enumerate(self.creators):
            if creator != other.creators[k]:
                return False

        # checking equality of created_date
        if self.created_date != other.created_date:
            return False

        # checking equality of modified_dates
        if len(self.modified_dates) != len(other.modified_dates):
            return False
        for k, modified_date in enumerate(self.modified_dates):
            if modified_date != other.modified_dates[k]:
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
            "created_date": self.created_date and self.created_date.isoformat(),
            "modified_dates": [
                mod_date.isoformat() for mod_date in self._modified_dates
            ],
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


class Creator:
    """Metadata for person who created an object.

    The creator has optional name, email and organisation properties.
    Separate given and family name values will be combined to a single
    value for name. (This prevents confusion due to different cultural
    conventions, see for example:
    https://uxmovement.com/forms/why-your-form-only-needs-one-name-field/)

    Parameters
    ----------
    name: str
        Full name of the creator. Optional, default None.
    email: str
        Email address of the creator. Optional, default None.
    organisation: str
        Name of the organisation of the creator, or the organisation that
        created the model. Optional, default None.
    """

    def __init__(
        self,
        name: Optional[str] = None,
        email: Optional[str] = None,
        organisation: Optional[str] = None,
        given_name: Optional[str] = None,
        family_name: Optional[str] = None,
    ):
        """Create an Creator metadata object."""
        self._name: Optional[str] = None
        self._email: Optional[str] = None
        self._organisation: Optional[str] = None

        self.name = Creator._fix_name(
            name=name, given_name=given_name, family_name=family_name
        )
        self.email = email
        self.organisation = organisation

    @staticmethod
    def _fix_name(
        name: Optional[str] = None,
        given_name: Optional[str] = None,
        family_name: Optional[str] = None,
    ):
        if name is not None and name != family_name:
            if given_name is not None or family_name is not None:
                raise ValueError(
                    """Too many name values were provided. Either a name or a
                    given and/or family name should be provided."""
                )
            return name
        if given_name is not None and family_name is not None:
            # This probably does not convert all names correctly.
            # Names should however preferentially be represented as a single value.
            return f"{given_name} {family_name}"
        elif given_name is not None:
            return given_name
        else:
            # This also covers the case where all values are None
            return family_name

    @property
    def name(self) -> Optional[str]:
        """Get the model creator name.

        Returns
        -------
        str
            Creator name.
        """
        return self._name

    @name.setter
    def name(self, value: Optional[str]) -> None:
        """Set the name of the model creator.

        Parameters
        ----------
        value: str
            Name of the creator.
        """
        self._name = value

    @property
    def email(self) -> Optional[str]:
        """Get the email address of the model creator.

        Returns
        -------
        str
            Email address.
        """
        return self._email

    @email.setter
    def email(self, value: Optional[str]) -> None:
        """Set the email of the model creator.

        Parameters
        ----------
        value: str
            Email address of the creator.
        """
        self._email = value

    @property
    def organisation(self) -> Optional[str]:
        """Get the organisation of the model creator.

        Returns
        -------
        str
            Organisation.
        """
        return self._organisation

    @organisation.setter
    def organisation(self, value: Optional[str]) -> None:
        """Set the organisation of the model creator.

        Parameters
        ----------
        value: str
            Organisation of the model creator.
        """
        self._organisation = value

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
            return Creator(**data)
        else:
            raise TypeError(f"Invalid format for Creator: {data}")

    def _asdict(self) -> Dict:
        d = {}
        if self.name is not None:
            d["name"] = self.name
        if self.email is not None:
            d["email"] = self.email
        if self.organisation is not None:
            d["organisation"] = self.organisation
        return d

    def to_dict(self) -> Dict:
        """Convert Creator to dictionary.

        Returns
        -------
        dict in this format
        {
            "name": str,
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
            f"('{self.name}', '{self.email}', "
            f"'{self.organisation}')"
        )

    def __eq__(self, other: Any) -> bool:
        """Determine whether the Creator is equal to another Creator object.

        Parameters
        ----------
        other: Creator
            Creator object to compare to.
        """
        if not isinstance(other, __class__):
            return False
        if self.name != other.name:
            return False
        if self.email != other.email:
            return False
        if self.organisation != other.organisation:
            return False
        return True
