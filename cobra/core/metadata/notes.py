# -*- coding: utf-8 -*-

from __future__ import absolute_import

import collections
import re

try:
    collectionsAbc = collections.abc
except AttributeError:
    collectionsAbc = collections


pattern_notes = re.compile(
    r"<(?P<prefix>(\w+:)?)p[^>]*>(?P<content>.*?)</(?P=prefix)p>",
    re.IGNORECASE | re.DOTALL
)


class Notes(collectionsAbc.MutableMapping):

    def __init__(self, notes_str=None):
        self._notes_str = notes_str
        self._data = {}
        self.update_notes_dict()

    def update_notes_dict(self):
        if self._notes_str and len(self._notes_str) > 0:
            notes_store = dict()
            for match in pattern_notes.finditer(self._notes_str):
                try:
                    # Python 2.7 does not allow keywords for split.
                    # Python 3 can have (":", maxsplit=1)
                    key, value = match.group("content").split(":", 1)
                except ValueError:
                    print("WARNING : Unexpected content format '{}'.",
                          match.group("content"))
                    continue
                self._data[key.strip()] = value.strip()
            return {k: v for k, v in notes_store.items() if len(v) > 0}
        else:
            return {}

    def update_notes_str(self, key, value):
        if self._notes_str is None:
            raise ValueError("Notes string is not a right place "
                             "to store key value pairs. Store them "
                             "at appropriate place in the document.")
        if not isinstance(value, str):
            print("WARNING : The value must be of type string. \n"
                  "Converting value in string type and "
                  "then putting in notes string.")
            value = str(value)
        pattern = re.compile(
            r"<(?P<prefix>(\w+:)?)p[^>]*>(\s*){}(\s*):(\s*)"
            r"(?P<content>.*?)(\s*)</(?P=prefix)p>".format(key),
            re.IGNORECASE | re.DOTALL
        )
        match = re.search(pattern, self._notes_str)
        if match is None:
            del self._data[key]
            raise ValueError("Notes string is not a right place to store "
                             "key value pairs. Store them at appropriate"
                             " place in the document.")
        else:
            start = match.start('content')
            end = match.end('content')
            modified_str = self._notes_str[:start] + \
                value + self._notes_str[end:]
        self._notes_str = modified_str

    def set_notes(self, notes_str):
        self._notes_str = notes_str
        self._data = {}
        if notes_str is None:
            return
        elif not isinstance(notes_str, str):
            raise TypeError("notes data must be of type "
                            "string: {}".format(notes_str))
        else:
            self.update_notes_dict()

    def get_notes_str(self):
        return self._notes_str

    def equals(self, new_notes):
        if self.get_notes_str() == new_notes.get_notes_str():
            return True
        else:
            return False

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
        return self._notes_str

    def __repr__(self):
        return self._notes_str
