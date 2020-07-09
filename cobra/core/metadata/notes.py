# -*- coding: utf-8 -*-

from __future__ import absolute_import

import collections
import re

try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections


class Notes(collectionsAbc.MutableMapping):
    """
    FIXME: documentation

    """
    # FIXME: should we abstract comments in the notes (?!)

    def __init__(self, notes_text: 'str' = None):
        self._notes_text = None
        self.notes_text = notes_text  # FIXME: rename notes_xhtml

    pattern_notes = re.compile(
        r"<(?P<prefix>(\w+:)?)p[^>]*>(?P<content>.*?)</(?P=prefix)p>",
        re.IGNORECASE | re.DOTALL
    )

    def update_notes_dict(self):
        """
        Updates notes dictionary according to key-value stored
        in notes string.
        """
        if self._notes_text and len(self._notes_text) > 0:
            for match in self.pattern_notes.finditer(self._notes_text):
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
        updated to store the new value.
        """
        # if notes string is empty
        if self._notes_text is None:
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
        match = re.search(pattern, self._notes_text)

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
            modified_str = self._notes_text[:start] + \
                value + self._notes_text[end:]
            self._notes_text = modified_str

    @property
    def notes_text(self):
        return self._notes_text

    @notes_text.setter
    def notes_text(self, value):
        if value is None:
            self._notes_text = value
            self._data = {}
        elif isinstance(value, str):
            self._notes_text = value
            self._data = {}
            self.update_notes_dict()
        else:
            raise TypeError("notes data must be of type "
                            "string: {}".format(value))

    def __eq__(self, other):
        if not isinstance(other, Notes):
            return False
        return self.notes_text == other.notes_text

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value
        self.update_notes_str(key, value)

    def __delitem__(self, key):
        del self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __str__(self):
        return self.notes_text

    def __repr__(self):
        return self.notes_text
