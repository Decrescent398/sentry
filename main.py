import spiceypy as sp
import numpy as np
from scipy.integrate import solve_ivp

from src.config import PLANTERY_METAKERNEL_TXT
from src.queries.spice_kernels import spice_setup, load_spice_kernels

from src.orbit import orbital_mechanics 
from src.orbit.methods import fitting, covariance, propagation
from src.prepare import residuals

#spice_setup()
load_spice_kernels()
sp.kclear()
sp.furnsh(str(PLANTERY_METAKERNEL_TXT))

def main():

    e_baseline, W, observations, x, t_start, t_end, trajectory_solution, start_index  = orbital_mechanics.get_residual_data()
    
    tolerance_step = 1e-8
    tolerance_grad = 1e-8
    tolerance_rel_Q = 1e-8

    Q_old = np.inf

    while True:
        
        e_baseline, W, Q, residual_states = residuals.calculate_residuals(observations, start_index, trajectory_solution)
        
        B = covariance.get_design_matrix(x, e_baseline, observations, t_start, t_end)
        
        C = B.T @ W @ B
        rhs = B.T @ W @ e_baseline
        
        lambda_ = 1e-3
        while True:
            C_lm = C + lambda_ * np.eye(6)
            dx = -np.linalg.solve(C_lm, rhs).flatten()
            
            x_trial = x + dx
            trial_trajectory = propagation.trajectory_solver(x_trial, t_start, t_end)
            _, _, Q_trial, _ = residuals.calculate_residuals(observations, start_index, trial_trajectory)
            
            if Q_trial < Q:
                x = x_trial
                lambda_ *=0.3
                break
            
            lambda_ *= 10.0
        
        step_norm = np.linalg.norm(dx)
        grad_norm = np.linalg.norm(rhs)
        
        if np.isfinite(Q_old):
            rel_Q = abs(Q_old - Q) / max(Q_old, 1.0)
        else:
            rel_Q = np.inf
        
        if step_norm < tolerance_step and grad_norm < tolerance_grad and rel_Q < tolerance_rel_Q:
            break
            
        Q_old = Q
        
    sigma = np.linalg.inv(C)
    
    return (x, t_start, t_end, residual_states)