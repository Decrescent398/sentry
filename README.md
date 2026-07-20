https://github.com/user-attachments/assets/b53e6e90-62de-40b6-aac5-fe22c79c363f

To run this, clone the repository, run pip install -r requirements.txt, and run notebooks/tests.ipynb in the stable branch (If you are using VSCode install the jupyter extension to run it).

This cannot be run on streamlit or any other jupyter notebook hosting service due to resources and timeouts, you need to download and run it only, takes about 3-4mins to run on most consumer laptops.

Significance of the results in the evals section of orbits.ipynb: 
  1. Mean residual is 0.13, this means that the calculated positions line up with the celestial coordinates of Eros to within a fraction of a human hair's width held at arm's length.
  2. Calculated position and calculated velocity of Eros (x) is within 30km and 10km/s of JPL Horizons calculations respectively.
