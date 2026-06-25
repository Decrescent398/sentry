import numpy as np

from src.orbit.methods.fitting import ObservationalError

def calculate_residuals(obs_data, start_index, trajectory_solution):
    observations = [obs for obs in obs_data[start_index+2:] if obs['stn'] != '270'] # omits and garbage observations from stn 270 and obs 1, 2
    result = [ObservationalError(obs, trajectory_solution) for obs in observations]

    residual_states, weights = zip(*result)
    residual_states = np.array(residual_states, dtype=np.float64)
    weights         = np.array(weights, dtype=np.float64)

    e_baseline = residual_states.flatten().reshape(-1, 1) # 2*sum(obs_index:), 1
    e_t = e_baseline.T # 1, 2*sum(obs_index:)

    W = np.diag(weights.flatten()) # 2*sum(obs_index:), 2*sum(obs_index:)
    q = e_t @ W @ e_baseline #1, 1
    
    return (e_baseline, W, observations)