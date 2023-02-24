#!/bin/bash

read -p "Run disk evol fortran script? " runevol
if [ $runevol = "y" ]
then
read -p "Enter the end time in millions of years [0.1 - 10]: " endtime
python3 $(find .. -name disk_input.py) $endtime
echo "Starting disk evol fortran scipt..."
./$(find .. -name diskevolfast)
fi
python3 $(find .. -name ice_anim.py)