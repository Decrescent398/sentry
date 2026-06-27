import numpy as np

from src.orbit.methods import fitting, propagation
from src.prepare import timesettings, residuals

def get_design_matrix(x_current, e_baseline, observation_list, t_start, t_end):
    relative_scale = 1e-8
    minimum_step = np.array([1e-3, 1e-3, 1e-3, 1e-6, 1e-6, 1e-6])
    
    num_residuals = len(e_baseline)
    params = 6
    B = np.zeros((num_residuals, params))
    
    start_index = timesettings.nearest_mpc_time_index(t_start, t_end)

    for j in range(params):
        
        h = np.maximum(relative_scale * np.abs(x_current[j]), minimum_step[j])
        
        x_plus = x_current.copy()
        x_minus = x_current.copy()
        
        x_plus[j] += h
        x_minus[j] -= h
        
        trajectory_plus = propagation.trajectory_solver(x_plus, t_start, t_end)
        e_plus, _, _, _ = residuals.calculate_residuals(observation_list, start_index, trajectory_plus)
        
        trajectory_minus = propagation.trajectory_solver(x_minus, t_start, t_end)
        e_minus, _, _, _ = residuals.calculate_residuals(observation_list, start_index, trajectory_minus)
        
        B[:, j] = ((e_plus - e_minus) / (2*h)).flatten()
    
    return B