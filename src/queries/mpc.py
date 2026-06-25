import requests
import numpy as np
import spiceypy as sp

from src.utils.units_and_conversions import AU_to_km, day_to_sec

def fetch_orbital_elements():
    response = requests.get("https://data.minorplanetcenter.net/api/get-orb", json={"desig": "Eros"})

    if response.ok:
        orb = response.json()[0]['mpc_orb']
    else:
        print("Error: ", response.status_code, response.content)
        
    orbital_elements = orb[0]['CAR']['coefficient_values']
    #coefficient values x, y, z, vx, vy, vz - state vector
    x = np.array(orbital_elements) 

    x[:3] *= AU_to_km #convert to km
    x[3:] *= AU_to_km / day_to_sec #convert to km/s

    t0_mjd = orb[0]['epoch_data']['epoch'] #mjd timestamp of state vector
    
    return (x, t0_mjd)

def fetch_observation_data():
    response = requests.get("https://data.minorplanetcenter.net/api/get-obs", json={"desigs": ["Eros"], "output_format": ["ADES_DF"]})

    if response.ok:
        obs_data = response.json()[0]['ADES_DF']
    else:
        print("Error: ", response.status_code, response.content)

    for n, obs in enumerate(obs_data): #convert to ephemeris time
        obs['obstime'] = sp.str2et(obs['obstime'])
        obs_data[n] = obs
        
    return obs_data