from src.queries import mpc
from src.prepare import timesettings, residuals, propagation_state
from src.orbit.methods import propagation

def get_residual_data():
    x, t0_mjd = mpc.fetch_orbital_elements()
    obs_data = mpc.fetch_observation_data()

    t_epoch, t_start, t_end = timesettings.set_initial_time_state(t0_mjd)

    x = propagation_state.rotate_x(x, t_start, t_epoch)
    trajectory_solution = propagation.trajectory_solver(x, t_start, t_end)

    start_index = timesettings.nearest_mpc_time_index(t_start, t_end)
    t_start = obs_data[start_index+2]
    
    e_baseline, W, observations = residuals.calculate_residuals(obs_data, start_index, trajectory_solution)

    return (e_baseline, W, observations, x, t_start, t_end)