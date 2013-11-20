#!/usr/bin/env python
assert __name__ == "__main__"
from os.path import abspath, split, join
import sys
sys.path.insert(0, abspath(join(split(__file__)[0], "..")))
from cobra.test import test_all
test_all()
