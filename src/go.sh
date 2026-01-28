#!/bin/bash

LIB_PATH="/home/basanta/lib/hmsbeagle:/home/basanta/lib/ncl"

export LD_LIBRARY_PATH="$LIB_PATH:$LD_LIBRARY_PATH"

export BEAGLE_LIB_DIR="/home/basanta/lib/hmsbeagle"

./strom "$@"
