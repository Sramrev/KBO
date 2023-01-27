#!/bin/bash

# Define the file name
fileName="script.py"

# Check if the file exists
if [ -f $fileName ]; then
  for i in $(seq 1 81)
  do
    # Create a new file with a part_i.py extension
    cp $fileName "part_$i.py"
    # Change the first line of the new file
    tail -n +2 "part_$i.py" > temp
    echo "part=$i" > "part_$i.py"
    cat temp >> "part_$i.py"
    rm temp
    echo "File part_$i.py has been created with first line as part=$i"
  done
else
  echo "The file does not exist."
fi
