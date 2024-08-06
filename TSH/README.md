 This folder contains the scripts and files for conducting the calculations
 needed for the F-tracking project

# 1. Models

To plot the PES (diabatic and adiabatic) for all the models, run:

    python plot_models.py

This will generate folders  HAM_0,  HAM_1, etc. with the figures

# 2. Methodologies 

In this project, I only use the GFSH method defined in the `recipes` folder. 


# 3a. Run individual TSH calculations

To run the TSH calculation, use the `tsh.py` file with the suitable arguments.
Run it e.g. like:

    python tsh.py --model_indx=0 --method_indx=0 --dt=4.0


# 3b. Run all TSH calculation

Edit the `run_all.py` script as needed, then:

    python run_all.py


# 4. Plot the population dynamics for all methods and models

Use Jupyter notebook: `plot-all.ipynb`
