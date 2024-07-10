#!/bin/bash

deactivate \
    && rm -rf ./crosslink-by-sequence-env \
    && python3 -m venv crosslink-by-sequence-env \
    && source ./crosslink-by-sequence-env/bin/activate \
    && pip install -e .