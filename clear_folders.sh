#!/bin/bash

# Define the directories to be cleared
DIR1="./data"
DIR2="./frames"

# Check if DIR1 exists and is a directory
if [ -d "$DIR1" ]; then
    echo "Clearing directory: $DIR1"
    rm -rf "$DIR1"/*
else
    echo "Directory $DIR1 does not exist."
fi

# Check if DIR2 exists and is a directory
if [ -d "$DIR2" ]; then
    echo "Clearing directory: $DIR2"
    rm -rf "$DIR2"/*
else
    echo "Directory $DIR2 does not exist."
fi

echo "Directories cleared."
