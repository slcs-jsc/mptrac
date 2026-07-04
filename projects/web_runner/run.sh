#! /bin/bash

# Create Python environment...
python3 -m venv venv
source ./venv/bin/activate
pip install -r requirements.txt

# Run flask with the local ERA-Interim test dataset...
python3 app.py --test
