import numpy as np
import spiceypy as sp

from src.utils.units_and_conversions import arcsec_to_rad, c_km_s
from src.utils.math_fns import angle_diff

from src.orbit.utils import station_utils
from src.orbit.utils import VFCC_utils

def get_observer_pos_j2000(stn, t_obs, properties):
    r_ecef = station_utils.stn_to_ecef(stn, properties)
        
    earth, _ = sp.spkezr('EARTH', t_obs, 'J2000', 'NONE', 'SUN')
    rot = sp.pxform('ITRF93', 'J2000', t_obs)
    r_obs_j2000 = rot @ r_ecef
    
    return earth[:3] + r_obs_j2000

def astrometric_error(details, dec_obs, arcsec_to_rad):
    
    sigma_ra  = None
    sigma_dec = None
    
    def rms(sigma_ra, sigma_dec):
        
        return {
            "rmsra"  : 1 / sigma_ra ** 2,
            "rmsdec" : 1 / sigma_dec ** 2}
            
    def largeMPCObserverError(details):
        
        sigma_ra  = float(details["rmsra"])
        sigma_dec = float(details["rmsdec"])
        
        if sigma_ra <= 0.2 and sigma_dec <= 0.2:
            return False
        else:
            return True
            
    if details["rmsra"] is not None and details["rmsdec"] is not None:
        
        sigma_ra  = float(details["rmsra"])
        sigma_dec = float(details["rmsdec"])
            
        if largeMPCObserverError(details):
            
            VFCCUncertainties = VFCC_utils.loadVFCC(details, dec_obs)
            
            sigma_ra  = VFCCUncertainties["sigma_ra"]
            sigma_dec = VFCCUncertainties["sigma_dec"]
            
    else:
        
        VFCCUncertainties = VFCC_utils.loadVFCC(details, dec_obs)
            
        sigma_ra  = VFCCUncertainties["sigma_ra"]
        sigma_dec = VFCCUncertainties["sigma_dec"]
        
    sigma_ra_rad  = sigma_ra * arcsec_to_rad
    sigma_dec_rad = sigma_dec * arcsec_to_rad
        
    rms_values = rms(sigma_ra=sigma_ra_rad, sigma_dec=sigma_dec_rad)
    
    return [rms_values["rmsra"], rms_values["rmsdec"]]

def ObservationalError(obs, trajectory_solution):

    t_obs   = obs['obstime']
    ra_obs  = np.deg2rad(float(obs['ra']))
    dec_obs = np.deg2rad(float(obs['dec']))
    
    state_at_obs = trajectory_solution(t_obs)
    r_asteroid = state_at_obs[:3]
    
    stn    = obs['stn']
    astcat = obs['astcat']
    
    stn_properties = station_utils.get_stn_properties(stn)
    obs_pos = 0.0
    
    if stn_properties['sin'] == None and stn_properties['cos'] == None:
        earth, _ = sp.spkezr('EARTH', t_obs, 'J2000', 'NONE', 'SUN')
        obs_pos = earth[:3]
    else:
        obs_pos = get_observer_pos_j2000(stn, t_obs, stn_properties)
        
    lt_old = 0.0 
        
    while True:
        
        rho = r_asteroid - obs_pos
        lt = np.linalg.norm(rho)
        
        if abs(lt - lt_old) < 1e-9:
            break
        
        lt_old = lt
        t_emit = t_obs - (lt / c_km_s)
        r_asteroid = trajectory_solution(t_emit)[:3]
        
    rho = r_asteroid - obs_pos
    rho_hat = rho / np.linalg.norm(rho)
    
    ra_pred = np.arctan2(rho_hat[1], rho_hat[0])
    ra_pred = np.mod(ra_pred, 2*np.pi)
    dec_pred = np.arcsin(rho_hat[2])
    
    residual = [angle_diff(ra_obs, ra_pred), dec_obs - dec_pred]
    
    rmsra = obs['rmsra']
    rmsdec = obs['rmsdec']
    
    telescope_details = {"stn": stn, "astcat": astcat, "rmsra": rmsra, "rmsdec": rmsdec}
    
    weight = astrometric_error(telescope_details, dec_obs, arcsec_to_rad)

    return residual, weight