#!/usr/bin/env python
# coding: utf-8
#
# This script provides an independent calculation of forward gravity from the 
# All American Canal MODFLOW model, for comparison with Heavy output.
# (test\aac\Python-Fortran_comparison.xlsx)
#
import os
import flopy
import flopy.utils.binaryfile as bf
import pandas as pd
from grav_module import timestep_gravity

print("dir {}".format(os.getcwd()))
os.chdir(r'../test/aac')
try:
    os.remove('sim_py.grav')
except:
    pass
    
model_basename = 'model1'
m = flopy.modflow.Modflow.load(model_basename+ '.nam', verbose=False, forgive=True)
hds = bf.HeadFile(model_basename+'.hds', precision = 'single')

# head is a [0][nr][nc] array, incrementing with the origin at the top left, downward and too the right.
# e.g., head[0][99] is the bottom row.
sy = m.upw.sy
reference_head = hds.get_data(totim=1)

# Get coordinate data and station names from this file
gravity_data = pd.read_csv('g_locations.csv', infer_datetime_format=True)
with open(os.path.join(m.model_ws,'sim_py.grav'),'a') as ofp4:
    ofp4.write('Station,TOTIM,g\n')
for station in set(gravity_data.Station):
    g_dic = {}       
    coord_data = gravity_data[gravity_data.Station.str.match(station)]
    coord_data = coord_data.reset_index(drop=True)
    x = coord_data['UTM_east'][0]
    y = coord_data['UTM_north'][0]
    z = coord_data['Elev'][0]

    ts = hds.get_times()
    for j in ts:
        print('Calculating gravity, observation: ' + station + ', time step: ' + str(j))
        timestep_head = hds.get_data(totim=int(j))
        timestep_g = timestep_gravity(m.modelgrid, sy.array, reference_head, timestep_head, x, y, z)
        g_dic.update({int(j):timestep_g})

    with open(os.path.join(m.model_ws,'sim_py.grav'),'a') as ofp4:
        for k, v in g_dic.items():
            ofp4.write('{},{},{:.2f}\n'.format(station, k, v))

print('grav target output file written')
