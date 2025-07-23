#!/bin/bash

# Check if directory argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

dir="$1"

# Check if the directory exists
if [ ! -d "$dir" ]; then
    echo "Directory not found: $dir"
    exit 1
fi

# Rename files ending with .n (n is one or more digits)
shopt -s nullglob
for file in "$dir"/*.[0-9]*; do
    [ -e "$file" ] || continue
    newname="${file%.[0-9]*}"
    mv "$file" "$newname"
done