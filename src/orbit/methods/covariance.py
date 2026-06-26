import numpy as np

from src.orbit.methods import fitting, propagation

def get_design_matrix(x_current, e_baseline, observation_list, t_start, t_end):
    num_residuals = len(e_baseline)
    params = 6
    B = np.zeros((num_residuals, params))

    perturbations = [1.0, 1.0, 1.0, 1e-4, 1e-4, 1e-4]

    for j in range(params):
        x_perturbed = x_current.copy()
        x_perturbed[j] += perturbations[j]
        
        perturbed_trajectory = propagation.trajectory_solver(x_perturbed, t_start, t_end)
        
        perturbed_residuals_list = []
        for obs in observation_list:
            res, _ = fitting.ObservationalError(obs, trajectory_solution=perturbed_trajectory)
            perturbed_residuals_list.extend(res)
            
        e_perturbed = np.array(perturbed_residuals_list, dtype=np.float64).reshape(-1, 1)
        
        B[:, j] = ((e_perturbed - e_baseline) / perturbations[j]).flatten()
    
    return B