# README
This is the code and data repository for the paper "Explicit Landau modeling of structural transitions in A2BB'O6 double perovskites reduces bias in critical temperature estimation." A breakdown of its contents is presented below:

## exp_data
Contains experimental x-ray diffraction data obtained from various sources, which are cited in the paper and included in "sources.txt". It contains 3 folders: one for raw data, another with unmerged crystallographic observations, and another with merged data, which was used for model fitting. The data merging script, processing.py, is also in this folder.

## qe
Contains the structures considered in DFT simulations and density of states results. The full output is too heavy and not the focus of the paper, so no other files are included

## results
Contains the coefficients obtained from the model fit: The folder for each material has one file that contains python output for the last Gauss-Newton step on both full-data (uncertain) and reduced-bias (certain) estimations. The results from full-data refinement are also aggregated here, in merged_results.txt and uncertainties.txt. The remaining file, T_shift.txt compares the transition temperature estimations (obtained using the coefficients here as input for arbitrary.py) with the indexing-system based estimations.

## scripts
 - **arbitrary.py**
 Calculates the behavior of q+, q-, epsilon for a given set of model parameters (vec, line 8) in the order  a+, T+, b+, a-, T-, b-, c, k, lambda+, lambda-, d+, d-, 0 (the last entry must be set to 0). The output is a file with a name defined in line 65, which has columns:  q+, q-, T, epsilon, 0,  status (0 for succesful convergence), and free energy.
 
 - **covars.py**
Analyzes the behavior of selected coefficients as a function of the Goldschmidt tolerance indices, and generates the temperature estimation comparison plots.

 -  **dos-qe.py**
Plots the electronic contribution to the internal energy as a function of temperature based on the results at qe/dos

- **estimate.py**
Performs the Gauss-Newton minimiaztion to find suitable coefficients for each material (chosen at line 6). The initial conditions are either defined by the user (line 44) or recovered from previous runs (commented-out lines 50-53) **Note: each step calls pyplot to generate a figure that gets  saved in the results folder and updated at each step, but the figures are not closed. Certain IDEs might open a new window on each step.** 
After num_iters steps have been performed (line 67), the results are stored in a npz file (line 147)

- **extrapolate.py**
Performs multiclass (lines 77-121) or binary (commented-out lines 169-190) classification trials, one algorithm at a time (lines 89-93/170-174) 

- **extrapolation.txt**
Contains the extrapolation dataset. Each line has entries: A, B, B', critical temperature, order reported in the literature, and doi of the source.

-**functions.py**
Contains definitions of free energy, its derivatives, gradient descent and Gauss-Newton iterations, which are used by arbitrary.py and estimate.py

-**ionradii.txt**
Contains ionic radii of elements as reported by [Shannon (1969)](https://doi.org/10.1107/S0567739476001551 "Open URL link"),  which are used to calculate tolerance factors.

-**plots.py**
Plots the results obtained by arbitrary.py.
