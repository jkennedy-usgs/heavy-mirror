Model files from:

Wildermuth, L., and J.R. Kennedy. 2022. MODFLOW-NWT
groundwater model demonstrating groundwater model calibration
with repeat microgravity measurements: U.S. Geological Survey 
data release. https://doi.org/10.5066/P9575C61

Gravity calculated at observation locations with the command:
$>Heavy.exe model1_hvy.nam

This writes the output files:
model1.lst
model1.out

Reference versions of these files are in the reference subdirectory.

Observation locations are defined in model1.hvy. Coordinates are
UTM Zone 11 N (NAD 83) relative to the model lower-left origin at
xll:658103
yll:3601515

The python file model1_run.py independently calculates forward 
gravity. Python output is compared to Heavy output in 
aac\Python-Fortran_comparison.xlsx.
