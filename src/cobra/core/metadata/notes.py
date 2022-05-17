import collections
import re
from typing import Iterator, Dict
from warnings import warn


class Notes(collections.MutableMapping):
    """Class representation of 'notes' of an object.

    The previous version of COBRApy was parsing entries of
    the form '<p> key : value </p>' and making
    a dict out of it. All other information was simply
    left out. When writing the model back to SBML, this
    dict was converted in the string:

        "<html xmlns = "http://www.w3.org/1999/xhtml" >
            <p>key1: value1</p>"
            <p>key2: value2</p>"
            ...
        </html>"

    The 'notes' attribute on an object stored this key: value
    dictionary.

    The current version of 'notes' has a dedicated class
    that behaves like a dict storing the key-value pairs
    present inside the notes string (making it backward
    compatible). In addition the the complete notes
    information is stored.

    The dict and the string of 'notes' are both synchronized
    with each other.
    Importantly, the 'notes' attribute is not meant to to
    store any machine-readable information. To enforce this
    behavior the addition of new key-values inside the
    'notes' dict is not permitted. Trying to do so
    will throw an ValueError. The KeyValuePairs should be used
    to store key:value information for objects.

    The complete 'notes' string is
    directly written to formats like "JSON", "YAML" etc when
    COBRA model is written in these format. And when writing
    SBML, 'notes' is initialized using the method:

            SBase.getNotesString()

    which makes the xhtml content of the notes using the string.

    Parameters
    ----------
    notes_xhtml : string
        The complete notes (xhtml) data in the form of a string.
    """

    # pattern checking for "<p> key : value </p>" type string
    PATTERN_PTAG = re.compile(
        r"<(?P<prefix>(\w+:)?)p[^>]*>(?P<content>.*?)</(?P=prefix)p>",
        re.IGNORECASE | re.DOTALL,
    )

    def __init__(self, notes_xhtml: str = None):
        self._data = {}
        self._notes_xhtml = None
        self.notes_xhtml = notes_xhtml

    @property
    def notes_xhtml(self) -> str:
        """Return the html content of notes in the form of a string."""
        return self._notes_xhtml

    @notes_xhtml.setter
    def notes_xhtml(self, value: str) -> None:
        """Set the notes_xhtml string"""
        if value is None:
            self._notes_xhtml = value
            self._data = {}
        elif isinstance(value, str):
            self._notes_xhtml = value
            self._data = {}
            self.update_notes_dict()
        else:
            raise TypeError(f"notes data must be of type string: {value}")

    def __eq__(self, other: "Notes") -> bool:
        if not isinstance(other, Notes):
            return False
        return self._notes_xhtml == other._notes_xhtml

    def __getitem__(self, key: str) -> str:
        return self._data[key]

    def __setitem__(self, key: str, value: str) -> None:
        if key not in self._data:
            raise ValueError(
                "Notes string is not a right place "
                "to store key value pairs. Store them "
                "at appropriate place in the document."
            )
        else:
            self._data[key] = value
            self.update_notes_str(key, value)

    def __delitem__(self, key: str) -> None:
        del self._data[key]

    def __iter__(self) -> Iterator:
        return iter(self._data)

    def __len__(self) -> int:
        return len(self._data)

    def __str__(self) -> str:
        if self._notes_xhtml is None:
            return ""
        return self.notes_xhtml

    def __repr__(self) -> str:
        return self.__str__()

    def update_notes_dict(self) -> None:
        """Updates notes dictionary according to key-value stored
        in notes string.
        """
        if self._notes_xhtml:
            for match in Notes.PATTERN_PTAG.finditer(self._notes_xhtml):
                try:
                    key, value = match.group("content").split(":", 1)
                except ValueError:
                    continue
                self._data[key.strip()] = value.strip()

    def update_notes_str(self, key: str, value: str) -> None:
        """Updates the notes string according to key-value pairs
        passed. If any such 'key' is present inside notes string
        having format '<p> key : oldvalue </p>', then it will be
        updated to store the new value. But if that 'key' is not
        present, an ValueError will be thrown.
        """
        # if notes string is empty
        if self._notes_xhtml is None:
            raise ValueError(
                "Notes string is not a right place "
                "to store key value pairs. Store them "
                "at appropriate place in the document."
            )

        # if value passed is not of type 'str'
        if not isinstance(value, str):
            warn(
                "The value must be of type string. \n"
                "Converting value to 'string' type and "
                "then putting in notes string...."
            )
            value = str(value)

        # pattern to search for inside notes string
        pattern = re.compile(
            rf"<(?P<prefix>(\w+:)?)p[^>]*>(\s*)"
            rf"{key}(\s*):(\s*)(?P<content>.*?)(\s*)</(?P=prefix)p>",
            re.IGNORECASE | re.DOTALL,
        )
        match = re.search(pattern, self._notes_xhtml)

        # if no such key-value substring is
        # already present inside notes string
        if match is None:
            del self._data[key]
            raise ValueError(
                "Notes string is not the right place "
                "to store key value pairs. Store them "
                "at appropriate place in the document."
            )
        # otherwise update the content
        else:
            start = match.start("content")
            end = match.end("content")
            modified_str = self._notes_xhtml[:start] + value + self._notes_xhtml[end:]
            self._notes_xhtml = modified_str

    @classmethod
    def notes_from_dict(cls, data_dict: Dict) -> "Notes":
        """Creates a new note based on a dictionary.

        This method can be used to create a completely new Notes object from a
        dictionary. It should be used when creating notes from scratch (such as import
        if the function already sets up a dict).

        This method will warn if using terms that should go into annotations.

        Parameters
        ----------
        data_dict: dict
            A dictionary that will be transformed to string.

        Returns
        -------
        Notes
            A new Notes object.
        """
        str_list = ['<html xmlns = "http://www.w3.org/1999/xhtml">']
        str_suffix = "</html>"
        for k, v in data_dict.items():
            str_list.append(f"<p>{k}: {v}</p>")
        str_list.append(str_suffix)
        return cls("\n".join(str_list))
