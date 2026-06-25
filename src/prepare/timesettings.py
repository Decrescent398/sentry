import spiceypy as sp
from datetime import datetime, timedelta

from src.queries import mpc

def set_initial_time_state(t0_mjd):
    jd_today = datetime.today().date() #todays date in jd

    t_epoch = sp.unitim(t0_mjd + 2400000.5, 'JDTDT', 'ET')  # convert mjd to jd
    t_start = sp.str2et(str(jd_today - timedelta(days=90))) # Last 90 days
    t_end = sp.str2et(str(jd_today))
    
    return (t_epoch, t_start, t_end)
    
def nearest_mpc_time_index(t_start, t_end):
    
    obs_data = mpc.fetch_observation_data()

    for n, obs in enumerate(obs_data):
        if obs['obstime'] >= t_start:
            t_start = obs_data[n]['obstime'] #propagating from nearest t_start datapoint
            t_end = obs_data[-1]['obstime'] #naturally t_end would be the last datapoint
            break
        
    total = len(obs_data)
    start_index = n
    
    return start_index