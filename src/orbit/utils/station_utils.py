import csv
import numpy as np

from src.config import STATIONS_CSV

from src.utils.units_and_conversions import R_earth

# stations.csv downloaded from https://www.minorplanetcenter.net/iau/lists/ObsCodes.html

def get_stn_properties(stn):
    
    with open(STATIONS_CSV, 'r', newline='\n') as stns:
        reader = csv.DictReader(stns)
        for row in reader:
            if row['Code'] == stn:
                return {
                    'Long.': float(row['Long.']) if row['Long.'].strip() != '' else None,
                    'cos': float(row['cos']) if row['cos'].strip() != '' else None,
                    'sin': float(row['sin']) if row['sin'].strip() != '' else None,
                }
    return None

def stn_to_ecef(stn, properties):
    
    lon = np.deg2rad(properties['Long.'])
    x = R_earth * properties['cos'] * np.cos(lon)
    y = R_earth * properties['cos'] * np.sin(lon)
    z = R_earth * properties['sin']
    
    return np.array([x, y, z])