# scripts/run_simulation.sh
#!/bin/bash

INPUT_FILE=$1
OUTPUT_DIR=$2

mkdir -p $OUTPUT_DIR

# Run Geant4 simulation
sim $INPUT_FILE
