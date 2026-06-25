import numpy as np
import spiceypy as sp
from scipy.integrate import solve_ivp

from src.orbit.methods.propagation import n_body_ode

def rotate_x(x, t_start, t_epoch):
    ecl_to_j2000 = sp.pxform('ECLIPJ2000', 'J2000', t_epoch)  #calculate state at t_start by propagating it from t_epoch to t_start
    x = np.concatenate([ecl_to_j2000 @ x[:3], ecl_to_j2000 @ x[3:]])
    x = solve_ivp(n_body_ode,
                    (t_epoch, t_start),
                    x,
                    method = "DOP853",
                    rtol=1e-12,
                    atol=1e-12,
                    dense_output=True,).y[:, -1]
    return x