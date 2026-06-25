import requests
import numpy as np
import spiceypy as sp

from src.utils.units_and_conversions import arcsec_to_rad
from main import main

x, t_start, t_end, residual_states = main()

print(np.sqrt(np.mean(residual_states**2)) / arcsec_to_rad)
print(np.percentile(np.abs(residual_states) / arcsec_to_rad, [50, 90, 95, 99]))

mean_residual = np.mean(np.abs(residual_states)) / arcsec_to_rad
print(f"Mean Residual in arcseconds: {mean_residual}")

command = "&COMMAND='433'"
obj_data = "&OBJ_DATA='YES'"
make_ephem = "&MAKE_EPHEM='YES'"
ephem_type = "&EPHEM_TYPE='VECTORS'"
center = "&CENTER='@10'" # Solar System Barycenter
ref_plane = "&REF_PLANE='F'" # BODY EQUATOR - Equatorial
out_units = "&OUT_UNITS='KM-S'"
start_time = f"&START_TIME='{sp.et2datetime(t_start).strftime("%Y-%m-%d %H:%M")}'"
stop_time = f"&STOP_TIME='{(sp.et2datetime(t_start) + timedelta(days=1)).strftime("%Y-%m-%d %H:%M")}'"
step_size = "&STEP_SIZE='1d'"

response = requests.get(f"https://ssd.jpl.nasa.gov/api/horizons.api?format=text{command}{obj_data}{make_ephem}{ephem_type}{center}{ref_plane}{out_units}{start_time}{stop_time}{step_size}")