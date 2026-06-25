VFCCLookupDefault = {
# (stn, astcat) : (ra, dec) arcsec 

# USNOA2.0 
('704', 'USNOA2')  : (0.63, 0.60),
('699', 'USNOA2')  : (0.62, 0.53),
('691', 'USNOA2')  : (0.30, 0.30),
('608', 'USNOA2')  : (0.61, 0.75),
('703', 'USNOA2')  : (0.69, 0.63),
('644', 'USNOA2')  : (0.29, 0.30),
('291', 'USNOA2')  : (0.46, 0.32),
('599', 'USNOA2')  : (0.39, 0.34),
('333', 'USNOA2')  : (0.55, 0.53),
('D35', 'USNOA2')  : (0.39, 0.38),

# USNOA1.0 
('704', 'USNOA1')  : (0.76, 0.73),
('691', 'USNOA1')  : (0.49, 0.46),

# USNOB1.0 
('699', 'USNOB1')  : (0.61, 0.54),
('644', 'USNOB1')  : (0.24, 0.20),
('691', 'USNOB1')  : (0.30, 0.28),
('291', 'USNOB1')  : (0.39, 0.26),

# UCAC1
('703', 'UCAC1')   : (0.63, 0.59),
('G96', 'UCAC1')   : (0.32, 0.27),
('E12', 'UCAC1')   : (0.50, 0.45),
('683', 'UCAC1')   : (0.79, 0.90),
('J75', 'UCAC1')   : (0.41, 0.37),
('106', 'UCAC1')   : (0.40, 0.39),
('143', 'UCAC1')   : (0.57, 0.47),

# UCAC2
('703', 'UCAC2')   : (0.63, 0.59),
('G96', 'UCAC2')   : (0.32, 0.27),
('E12', 'UCAC2')   : (0.50, 0.45),
('683', 'UCAC2')   : (0.79, 0.90),
('J75', 'UCAC2')   : (0.41, 0.37),
('106', 'UCAC2')   : (0.40, 0.39),
('143', 'UCAC2')   : (0.57, 0.47),

# Gaia2
('T14', 'Gaia2')   : (0.10, 0.10),
('T12', 'Gaia2')   : (0.10, 0.10),
('T09', 'Gaia2')   : (0.10, 0.10),
('Y28', 'Gaia2')   : (0.30, 0.30),
('568', 'Gaia2')   : (0.10, 0.10),
('G83', 'Gaia2')   : (0.20, 0.20),
('309', 'Gaia2')   : (0.20, 0.20),

# Gaia3
('T14', 'Gaia3')   : (0.10, 0.10),
('T12', 'Gaia3')   : (0.10, 0.10),
('T09', 'Gaia3')   : (0.10, 0.10),
('Y28', 'Gaia3')   : (0.30, 0.30),
('568', 'Gaia3')   : (0.10, 0.10),
('G83', 'Gaia3')   : (0.20, 0.20),
('309', 'Gaia3')   : (0.20, 0.20),

# Gaia3E
('T14', 'Gaia3E')  : (0.10, 0.10),
('T12', 'Gaia3E')  : (0.10, 0.10),
('T09', 'Gaia3E')  : (0.10, 0.10),
('Y28', 'Gaia3E')  : (0.30, 0.30),
('568', 'Gaia3E')  : (0.10, 0.10),
('G83', 'Gaia3E')  : (0.20, 0.20),
('309', 'Gaia3E')  : (0.20, 0.20),

# Tycho-2
('689', 'Tycho-2') : (0.20, 0.21),

}

VFCCAstcat = {
# (astcat) : (ra, dec) arcsec 

'Tycho-2'          : (0.24, 0.25),
'UCAC2'            : (0.53, 0.49),
'UCAC1'            : (0.53, 0.49),
'UCAC4'            : (0.30, 0.30),
'USNOB1'           : (0.48, 0.42),
'USNOA1'           : (0.72, 0.69),
'USNOA2'           : (0.61, 0.58),
'Gaia2'            : (0.20, 0.20),
'Gaia3'            : (0.20, 0.20),
'Gaia3E'           : (0.20, 0.20),
'ATLAS2'           : (0.20, 0.20),

}

VFCCStn = {
# (stn) : (ra, dec) arcsec 

'645'              : (0.30, 0.30),
'673'              : (0.30, 0.30),
'689'              : (0.50, 0.50),
'950'              : (0.50, 0.50),
'H01'              : (0.30, 0.30),
'J04'              : (0.40, 0.40),
'W84'              : (0.50, 0.50),
'LCO'              : (0.40, 0.40),

}

VFCCLookup = {
    
    "Default" : VFCCLookupDefault,
    "Astcat"  : VFCCAstcat,
    "Station" : VFCCStn,
    
}

def loadVFCC(details, dec_obs, lookup=VFCCLookup):
    
    if (details["stn"], details["astcat"]) in VFCCLookup["Default"]:
        
        sigma_ra  = VFCCLookup["Default"][(details["stn"], details["astcat"])][0]
        sigma_dec = VFCCLookup["Default"][(details["stn"], details["astcat"])][1]
        
    elif details["stn"] in VFCCLookup["Station"] and details["astcat"] in VFCCLookup["Astcat"]:
        
        sigma_ra  = VFCCLookup["Astcat"][details["astcat"]][0]
        sigma_dec = VFCCLookup["Astcat"][details["astcat"]][1]
        
        if sigma_ra < VFCCLookup["Station"][details["stn"]][0]:
            sigma_ra = VFCCLookup["Station"][details["stn"]][0]
            
        if sigma_dec < VFCCLookup["Station"][details["stn"]][1]:
            sigma_dec = VFCCLookup["Station"][details["stn"]][1]
        
    elif details["astcat"] in VFCCLookup["Astcat"]:
        
        sigma_ra  = VFCCLookup["Astcat"][details["astcat"]][0]
        sigma_dec = VFCCLookup["Astcat"][details["astcat"]][1]
        
    elif details["stn"] in VFCCLookup["Station"]:
        
        sigma_ra  = VFCCLookup["Station"][details["stn"]][0]
        sigma_dec = VFCCLookup["Station"][details["stn"]][1]
        
    else:
        
        sigma = 0.2
        
        sigma_ra  = sigma * np.cos(dec_obs)
        sigma_dec = sigma
        
    return {"sigma_ra": sigma_ra, "sigma_dec": sigma_dec}