#!/bin/bash
#
#   Run this script to test scROB.
#   
#   2023/05/04 @yanisaspic

Rscript src/tests/style.R
black .

mypy . --ignore-missing-imports
python3 -m unittest
