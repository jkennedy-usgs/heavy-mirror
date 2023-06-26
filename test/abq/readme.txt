Model files from:
Myers, N.C., and Friesz, P.J., 2019, Hydrogeologic framework and delineation of transient
areas contributing recharge and zones of contribution to selected wells in the upper 
Santa Fe Group aquifer, southeastern Albuquerque, New Mexico, 1900–2050: U.S. Geological 
Survey Scientific Investigations Report 2019–5052, 73 p., 
https://doi.org/10.3133/sir20195052.

Model downloaded from https://water.usgs.gov/GIS/metadata/usgswrd/XML/sir2019-5052.xml

Zone arrays are defined in the externalfiles/arrays directory.

Gravity calculated at a grid of observation locations with the command:
$>Heavy.exe -1 -g 50 model1_hvy.nam

This writes the output files:
tran_hvy.lst
tran_hvy.out

Reference versions of these files are in the reference subdirectory.

This example demonstrates calculating gravity change relative to the first time step
using the -t1 flag. In this model the initial head is not steady-state.