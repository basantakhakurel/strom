#!/bin/bash

# setup variables
BASE_DIR=$(pwd)
SRC_DIR="$BASE_DIR/src"
BUILD_DIR="$BASE_DIR/build"
INSTALL_DIR="$BASE_DIR/distr"

# Create the build dir if it does not already exist
if [ ! -d "$BUILD_DIR" ]; then
    mkdir "$BUILD_DIR"
fi

cd "$BASE_DIR"

meson setup "$BUILD_DIR" "$SRC_DIR" --prefix="$INSTALL_DIR"

ninja -C "$BUILD_DIR"

ninja -C "$BUILD_DIR" install

cd "$BASE_DIR"
