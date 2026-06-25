import spiceypy as sp
import numpy as np
from scipy.integrate import solve_ivp

from src.config import PLANTERY_METAKERNEL_TXT
from src.queries.spice_kernels import spice_setup, load_spice_kernels

from src.orbit import orbital_mechanics 
from src.orbit.methods import fitting, covariance, propagation

# spice_setup()
load_spice_kernels()
sp.kclear()
sp.furnsh(str(PLANTERY_METAKERNEL_TXT))

def main():

    e_baseline, W, observations, x, t_start, t_end  = orbital_mechanics.get_residual_data()

    max_iterations = 8
    tolerance_dx = 1e-3   
    tolerance_dQ = 1e-2   

    Q_old = float('inf')

    for iteration in range(max_iterations):
        
        trajectory_solution = propagation.trajectory_solver(x, t_start, t_end)
        
        result = [fitting.ObservationalError(obs, trajectory_solution=trajectory_solution) for obs in observations]
        residual_states, _ = zip(*result)
        residual_states = np.array(residual_states, dtype=np.float64)
        
        e_baseline = residual_states.flatten().reshape(-1, 1)
        
        Q = (e_baseline.T @ W @ e_baseline).item()
        
        B = covariance.get_design_matrix(x, e_baseline, observations, t_start, t_end)
        
        C = B.T @ W @ B
        rhs = B.T @ W @ e_baseline
        
        damping_factor = 1e-3
        C_stabilized = C + damping_factor * np.eye(6)
        
        dx = -np.linalg.solve(C_stabilized, rhs).flatten()
        
        x += dx
        
        pos_step_mag = np.linalg.norm(dx[:3])
        
        if pos_step_mag < tolerance_dx or abs(Q_old - Q) < tolerance_dQ:
            break
            
        Q_old = Q
        
    sigma = np.linalg.inv(C)
    
    return (x, t_start, t_end, residual_states)