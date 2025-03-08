#!/bin/bash

DYLD_LIBRARY_PATH="/usr/local/lib:$DYLD_LIBRARY_PATH" ./strom $@
