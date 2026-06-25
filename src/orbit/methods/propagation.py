import numpy as np
import spiceypy as sp
from scipy.integrate import solve_ivp

from src.utils.units_and_conversions import sun_mu

def n_body_ode(t, state):
    
    eros_position = state[:3]
    eros_velocity = state[3:]
    
    r = eros_position
    v = eros_velocity
    
    r_norm = np.linalg.norm(r)
    a = -sun_mu * r / r_norm**3
    
    ssbs = [
        # (name, mu)
        
        # Terrestrial
        ('1', 22031.78000000002), # MERCURY
        ('2',   324858.592), # VENUS
        ('3',   398600.435), # EARTH
        ('4',    42828.375), # MARS
        
        # Satellites
        ('301',               4902.800), # MOON
        
        # Asteroids
        ('2000001',              63.129999999999995),
        ('2000002',              13.73),
        ('2000004',              17.28999999999999),
        
        ('2000010',            5.78),               # Hygiea
        ('2000015',            2.10),               # Eunomia
        ('2000016',            1.81),               # Psyche
        ('2000029',            0.86),               # Amphitrite
        ('2000052',            1.5899999999999999), # Europa
        ('2000065',            0.9099999999999999), # Cybele
        ('2000087',            0.9899999999999999), # Sylvia
        ('2000088',            1.02),               # Thisbe
        ('2000511',            2.259999999999999),  # Davida
        ('2000704',            2.189999999999999),  # Interamnia
        
        # Jovian
        ('5', 126712764.1), # JUPITER
        ('6',  37940585.2), # SATURN
        ('7',  5794548.600000008), # URANUS
        ('8', 6836527.100580023), # NEPTUNE
    ]
    
    for name, mu in ssbs:
        
        # i is planet, j is eros
        
        planet, _ = sp.spkezr(name, t, 'J2000', 'NONE', 'SUN') #calculates state vector of planet
        
        r_i = planet[:3]
        r_i_norm = np.linalg.norm(r_i)
        
        r_ji = r - r_i
        r_ji_norm = np.linalg.norm(r_ji)
        
        # subtracting pull due to sun
        a += -mu * (r_ji / r_ji_norm**3 + r_i / r_i_norm**3)
    
    return np.concatenate((v, a)) # vx, vy, vz, ax, ay, az

def trajectory_solver(x, t_start, t_end):

    trajectory_solution = solve_ivp(n_body_ode, #precompute entire trajectory from t_start to t_end
                                    (t_start, t_end),
                                    x,
                                    method = "DOP853",
                                    rtol=1e-12,
                                    atol=1e-12,
                                    dense_output=True,).sol
    
    return trajectory_solution