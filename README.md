# Dose_optimization

1. Load_raw_data.m – reads in the excel files from each experiment ( for all cell types and treatment regimens)- generates trajraw.mat
2. Filter_data_231.m loads in trajraw.mat and selects only the 231 cells from this structure, and then selects only those that have only received one treatment. Assigns a color and fits a bi-exponential model to each individual cell well. Records time to reach 2* baseline (critical time).  Generates trajfit231.mat
3. Fit_N_231.m loads in trajfit231.mat and combines wells from the same treatment scenario to get a mean and standard deviation vector for each treatment. Also truncates the data just below the carrying capacity, and adds an effective dose U(t) based on the concentration of dox for each treatment. Then fits the N_t data to the model using the fit_fxn_Greene.m which has 3 ways of calibrating- using normpdf fxn, using weighted least-squares, and using just least-squares. Generated pfitN and trajsumfit231.mat
4. Fit_N_phi_231.m loads in trajsumfit231.mat to get the N(t) data and phi_t_est_pyth.csv to get the phi(t) data. Fits the joint data using fit_fxn_N_phi.m which does a fit using both the weighted relative error in N and the weighted relative error in phi. It also does the fitting using just the N data and just the weighted relative error in N. Bother outputs of parameter estimates and models are saved.
5. Test_param_predictions loads in trajsumfit231 and the fit parameters from integrated fit and fit on N alone. It then runs everything forward with those parameters- comparing the fit to N(t) for the integrated fit and the N(t) alone fit, and comparing the predictions for each. 





Parameters for traj structure

Never changes
time				      raw values
rawN				      raw values
N0true				    rawN(1)
welllabel			    raw well string
well				      well string
column				    column # (STRING)
accdose				    dose + prevdose

Never changes (right now)
drug				      currently all dox
doseduration			currently all 24

Changes: applies to entire plate
celltype		      defaults to MCF-7
date				      date string
tdose				      aspiration time of dose
seed				      default 2000, one data set has 3000

Changes: may change between wells
dose				      dose amount
dosenum				    which # dose the current dose is
WPT				        weeks post treatment
prevdose			    all previous dosage amounts (may be vector or 0)
doseints			    all previous dosage intervals (may be vector or 0)
numdoses			    total # of doses including current

If dose = 0, WPT/dosenum/doseints/numdoses are still treated as if dosage was present
