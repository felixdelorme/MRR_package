# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 11:50:12 2022

@author: rebouda
"""

#Commons
import netCDF4 as nc
import pandas as pd
import numpy as np
import os

path='/Users/felixdelorme/ESISAR/Stage/stage_felix_2024_MRR_Bolivie/MK_packaged_v2/MK_processed/PYRAMIDE/test/'
file1='0508_60s_PLATAFORMA.nc'
filelist= os.listdir(path)

for file in filelist:
    ds=nc.Dataset(path+file, mode='r')
    time=ds.variables['time']
    print(os.path.basename(file))
    print('start:',pd.to_datetime(time[:], unit='s')[0], 
          '\nend:',pd.to_datetime(time[:], unit='s')[-1])
    ds.close()

"""mergefile='/Users/felixdelorme/ESISAR/Projets/Stage/stage_felix_2024_MRR_Bolivie/CLAUDIO.DURAN/MK_packaged_v2/Data/PLATAFORMA/MK_processed/202309/202309.nc'

ds=nc.Dataset(mergefile, mode='r')
time=ds.variables['time']
pd.to_datetime(time[1430:1450], unit='s')
print(os.path.basename(file))
print('start:',pd.to_datetime(time[:], unit='s')[0], 
      '\nend:',pd.to_datetime(time[:], unit='s')[-1])
ds.close()"""

