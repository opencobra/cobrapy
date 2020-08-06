# -*- coding: utf-8 -*-

import collections
import re

try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections


class Notes(collectionsAbc.MutableMapping):
    """ Class representation of 'notes' of COBRA object.

    The previous version of COBRApy was parsing out any-
    thing of the form '<p> key : value </p>' and making
    a dict out of it. All other information was simply
    left out. And while writing model back to SBML, this
    dict was converted in the string:

        "<html xmlns = "http://www.w3.org/1999/xhtml" >
            <p>key1: value1</p>"
            <p>key2: value2</p>"
            ...
        </html>"

    So in the name of 'notes', we were having a dict
    containing these key-value pairs in COBRApy in
    previous version.

    The current version of 'notes' have a dedicated class
    that behaves like a dict storing the key-value pairs
    present inside the notes string (making it backward
    compatible) and also stores the complete notes data
    in the form of a string inside it.

    The dict and the string of 'notes' are both synchronized
    with each other. When a notes string is set, the dict will
    get updated according to the key-value pairs present inside
    the notes string. If you change a value corresponding to a
    key inside the notes dict, the notes string will also be
    updated. However, one thing to note is, since the 'notes'
    attribute is not a right place to store any machine
    readable information and only human-readable information
    should be put into it, so addition of any new key-values
    inside the 'notes' dict is not allowed. Trying to do so
    will throw an ValueError. The complete 'notes' string is
    directly written to formats like "JSON", "YAML" etc when
    COBRA model is written in these format. And when writing
    SBML, 'notes' is initialized using the method:

            SBase.getNotesString()

    which makes the xhtml content of the notes using the string.

    Parameters
    ----------
    notes_xhtml : string
        The complete notes (xhtml) data in the form of a string.

    Attributes
    ----------
    notes_xhtml : string
        The string corresponding to the notes (xhtml) data

    """
    # FIXME: should we abstract comments in the notes (?!)

    def __init__(self, notes_xhtml: 'str' = None):
        self._notes_xhtml = None
        self._data = {}
        self.notes_xhtml = notes_xhtml

    pattern_notes = re.compile(
        r"<(?P<prefix>(\w+:)?)p[^>]*>(?P<content>.*?)</(?P=prefix)p>",
        re.IGNORECASE | re.DOTALL
    )

    def update_notes_dict(self):
        """
        Updates notes dictionary according to key-value stored
        in notes string.
        """
        if self._notes_xhtml and len(self._notes_xhtml) > 0:
            for match in self.pattern_notes.finditer(self._notes_xhtml):
                try:
                    # Python 2.7 does not allow keywords for split.
                    # Python 3 can have (":", maxsplit=1)
                    key, value = match.group("content").split(":", 1)
                except ValueError:
                    continue
                self._data[key.strip()] = value.strip()

    def update_notes_str(self, key, value):
        """
        Updates the notes string according to key-value pairs
        passed. If any such 'key' is present inside notes string
        having format '<p> key : oldvalue </p>', then it will be
        updated to store the new value. But if that 'key' is not
        present, an ValueError will be thrown.
        """
        # if notes string is empty
        if self._notes_xhtml is None:
            raise ValueError("Notes string is not a right place "
                             "to store key value pairs. Store them "
                             "at appropriate place in the document.")

        # if value passed is not of type 'str'
        if not isinstance(value, str):
            print("WARNING : The value must be of type string. \n"
                  "Converting value to 'string' type and "
                  "then putting in notes string....")
            value = str(value)

        # pattern to search for inside notes string.
        pattern = re.compile(
            r"<(?P<prefix>(\w+:)?)p[^>]*>(\s*){}(\s*):(\s*)"
            r"(?P<content>.*?)(\s*)</(?P=prefix)p>".format(key),
            re.IGNORECASE | re.DOTALL
        )
        match = re.search(pattern, self._notes_xhtml)

        # if no such key-value substring is
        # already present inside notes string
        if match is None:
            del self._data[key]
            raise ValueError("Notes string is not a right place "
                             "to store key value pairs. Store them "
                             "at appropriate place in the document.")
        # otherwise update the content
        else:
            start = match.start('content')
            end = match.end('content')
            modified_str = self._notes_xhtml[:start] + \
                value + self._notes_xhtml[end:]
            self._notes_xhtml = modified_str

    @property
    def notes_xhtml(self):
        """ Return the html content of notes in the form
        of a string.
        """
        return self._notes_xhtml

    @notes_xhtml.setter
    def notes_xhtml(self, value):
        """ Set the notes_xhtml string
        """
        if value is None:
            self._notes_xhtml = value
            self._data = {}
        elif isinstance(value, str):
            self._notes_xhtml = value
            self._data = {}
            self.update_notes_dict()
        else:
            raise TypeError("notes data must be of type "
                            "string: {}".format(value))

    def __eq__(self, other):
        if not isinstance(other, Notes):
            return False
        return self.notes_xhtml == other.notes_xhtml

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        if key not in self._data:
            raise ValueError("Notes string is not a right place "
                             "to store key value pairs. Store them "
                             "at appropriate place in the document.")
        else:
            self._data[key] = value
            self.update_notes_str(key, value)

    def __delitem__(self, key):
        del self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __str__(self):
        return self.notes_xhtml

    def __repr__(self):
        return self.notes_xhtml
