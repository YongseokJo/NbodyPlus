import pandas as pd
import numpy as np

column_name_nbody = ['M', 'x', 'y', 'z', 'vx', 'vy', 'vz'] # Msol, pc, km/s
df = pd.read_csv('./astro.NBODY', sep='\s+', names=column_name_nbody)
                 
column_name_enzo = ['x', 'y', 'z', 'vx', 'vy', 'vz', 'M'] # kpc, km/s, 1e9 Msol
enzo_file = pd.DataFrame(index=range(df.shape[0]), columns=range(df.shape[1]))
enzo_file.columns = column_name_enzo

enzo_file['x'] = df['x']*1e-3
enzo_file['y'] = df['y']*1e-3
enzo_file['z'] = df['z']*1e-3
enzo_file['vx'] = df['vx']
enzo_file['vy'] = df['vy']
enzo_file['vz'] = df['vz']
enzo_file['M'] = df['M']*1e-9
                 
enzo_file.to_csv('./nbody.dat', sep=' ', index=False, header=False)
