#!/bin/sh
module use /projects/academic/chrest/modules
module load chrest/release
./exe --useYaml=true --chemfile=../tchemExample/inputs/LL2KGB_AllRange_UnitsTesting.yaml --samplefile=../tchemExample/inputs/sample2.dat