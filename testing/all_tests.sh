#!/bin/sh
cd ABLATEtest/2S_CH4_CM2
rm _ignitionDelay2S_CH4_CM2/ignitionDelayTemperature.txt
echo 'testing ABLATE with 2S_CH4_CM2'
bash 2S_CH4_CM2.sh
echo ''

cd ../pentane
echo 'testing ABLATE with pentane'
bash pentane.sh
echo ''

cd ../../DAKOTAtest/driver_success
rm _ignitionDelay2S_CH4_CM2/ignitionDelayTemperature.txt
echo 'testing a block of code from ABLATE-DAKOTA driver with 2S_CH4_CM2'
bash driver_success.sh
echo ''

cd ../driver_fail
echo 'testing full ABLATE-DAKOTA driver with 2S_CH4_CM2'
bash driver_fail.sh
echo ''

cd ../../TCHEMtest/tchemExampleBuild
echo 'testing TCHEM with 32 samples of pentane'
bash pentane_32_samples.sh
echo ''

echo 'testing TCHEM with 1 sample of pentane'
bash pentane_1_sample.sh
echo '\n'

cd ../../../
git status