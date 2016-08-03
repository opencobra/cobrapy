#!/bin/bash

pip install pip --upgrade
pip install numpy scipy python-libsbml cython coveralls jsonschema six matplotlib pandas
wget https://opencobra.github.io/pypi_cobrapy_travis/esolver.gz
gzip -d esolver.gz
chmod +x esolver
export PATH=$PATH:$PWD
mkdir -p $HOME/.config/matplotlib
echo 'backend: Agg' >> $HOME/.config/matplotlib/matplotlibrc
python setup.py develop
coverage run --source=cobra setup.py test
