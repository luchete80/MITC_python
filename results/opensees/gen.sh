#!/bin/bash

# Root folder of OpenSees source
SRC_ROOT="/home/weldform-pc/Numerico/OpenSees/SRC"

echo "include_directories("

# Find all unique directories containing header files
find "$SRC_ROOT" -type d | sort | while read -r dir; do
    echo "    \${PROJECT_SOURCE_DIR}/$dir"
done

echo ")"
