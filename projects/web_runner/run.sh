#! /bin/bash

# Create Python environment...
python3 -m venv venv
source ./venv/bin/activate
pip install -r requirements.txt

# Run flask...
python3 app.py
