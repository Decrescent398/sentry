To run this, clone the repository, run pip install -r requirements.txt, and run notebooks/orbits.ipynb (If you are using VSCode install the jupyter extension to run).

This cannot be run on a streamlit or any other nb hosting service, you need to download and run it, takes about 3-4mins to run

Significance of the results in the evals section of orbits.ipynb: 
  1. Mean residual is 0.13, this means that the calculated positions line up with the celestial coordinates of Eros to within a fraction of a human hair's width held at arm's length.
  2. Calculated position of Eros (x) is within 10s of kms of JPL Horizons calculations as displayed.
