import os
from datetime import datetime

from src.config import KERNEL_DIR, EXAMPLE_SPICEKERNEL, EARTH_LATEST_HIGH_PRESCISION

def spice_setup():
    os.remove(EXAMPLE_SPICEKERNEL)
    
    today = datetime.today().strftime("%D")
    day = today.split('/')[1]
    reset = day % 30
    
    if reset: 
        os.remove(EARTH_LATEST_HIGH_PRESCISION)
    return 
        

def load_spice_kernels():
    
    # The meta kernel file contains entries pointing to the following SPICE kernels, which the user needs to download.
    #   https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls
    #   https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp
    #   https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc
    #   https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc
    #   https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/codes_300ast_20100725.tf
    #   https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc
    #   https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/codes_300ast_20100725.bpc

    #   The following is the contents of a metakernel that was saved with
    #   the name 'planetaryMetaK.txt'.
    #   \begindata
    #   KERNELS_TO_LOAD=(
    #     'kernels/naif0012.tls',
    #     'kernels/de440.bsp',
    #     'kernels/earth_latest_high_prec.bpc',
    #     'kernels/gm_de431.tpc',
    #     'kernels/codes_300ast_20100725.tf',
    #     'kernels/pck00010.tpc',
    #     'kernels/codes_300ast_20100725.bpc'
    #   )
    #   \begintext
    
    kernels = {
    "lsk/": ['naif0012.tls',],
    
    "spk/": ['planets/de440.bsp', 
            'asteroids/codes_300ast_20100725.tf',
            'asteroids/codes_300ast_20100725.bsp',],
    
    "pck/": ['earth_latest_high_prec.bpc', 
            'gm_de431.tpc',
            'pck00010.tpc',]
    }
    
    url = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/'
    kernels_loaded = os.listdir(KERNEL_DIR)

    for kernel_type, items in kernels.items():
        for item in items:
            item_name = item.split('/')[-1]
            if item_name in kernels_loaded:
                continue
            else:
                download_url = url + kernel_type + item
                urllib.request.urlretrieve(download_url, f"{KERNEL_DIR}/{item_name}")      
    return