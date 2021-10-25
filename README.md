# Cam_Clay_Explicit_Undrained

MCC_Explicit_mono code is to run an undrained monotonic test for given preconsolidation pressure and current mean confining stress.
Typical input lines in the code look like:

This program is for predicting soil's behavior using MCC model in a Undrained TRIAXIAL TEST\\
Enter the maximum Isotropic consolidation stress (Po')(kPa):392
Please define the initial stress state in p-q space-
Enter the initial mean effective confining stress on the sample (kPa):305
Running Triaxial Undrained Test
Please Enter the strain level you wish to plot the results(%):8
Enter the Strain Increament(%):0.01

MCC_Explicit_cyc code is to run an cyclic undrained test for given preconsolidation pressure and current mean confining stress.
Typical input lines in the code look like:

This program is for predicting soil's behavior using MCC model in a Cyclic Undrained TRIAXIAL TEST
Enter the maximum Isotropic consolidation stress (Po')(kPa):392
Please define the initial stress state in p-q space-
Enter the initial mean effective confining stress on the sample (kPa):392
Enter peak-peak deviotric stress for cyclic loading:200
Number of cyclic loadings:5
Running Cyclic Triaxial Undrained Test
Enter the Strain Increament(%):0.01

