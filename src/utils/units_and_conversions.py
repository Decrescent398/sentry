import numpy as np

AU_to_km = 1.495978707e8 # AU to Km
day_to_sec = 86400.0 # Day to Seconds
arcsec_to_rad = np.deg2rad(1/3600) # Arcsecond to Radian

sun_radius = 695700 # Radius of the Sun (km)
sun_mu = 1.32712440018e11  # Gravitational parameter of the Sun (km^3/s^2)

R_earth = 6378.137 # Radius of the Earth(km)
R_earth_polar = 6356.7523 # Polar radius of the Earth(km)

c_km_s = 299792.458 #Speed of light (s)