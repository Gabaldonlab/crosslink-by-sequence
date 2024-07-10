#!/bin/bash
deactivate && rm -rf testy && python3 -m venv testy && source testy/bin/activate && pip install biopython pandas && python3 -m crosslink_by_sequence