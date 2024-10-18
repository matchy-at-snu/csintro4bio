#!/bin/bash -e

# get the location of this file
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Run all the scripts in the correct order
for i in $(seq 0 1 5); do
    python "$DIR/mission/0$i/Mission$i.py"
done
