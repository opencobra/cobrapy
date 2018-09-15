# -*- coding: utf-8 -*-

"""Provide JSON schemata that can be used for validation."""

from __future__ import absolute_import

import json

from importlib_resources import open_text


with open_text("cobra.io.schemata", "model_schema.json",
               encoding="utf-8") as file_handle:
    MODEL_SCHEMA = json.load(file_handle)
