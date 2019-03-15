# Dose_optimization

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
