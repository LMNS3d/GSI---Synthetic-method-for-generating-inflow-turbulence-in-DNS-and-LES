# GSI: A generalized synthetic method for generating inhomogeneous inflow turbulence in DNS and LES

This repository covers the MATLAB implementation of the GSI for generating the inflow boundary condition
for large-eddy simulation or direct numerical simulation of the homogeneous and inhomogeneous turbulence,
which is proposed in Ref. (Ying, Li & Fu, Generalized synthetic inflow generation method for divergence-free 
inhomogeneous turbulence, JCP, 2025).

The examples correspond to the simulated cases presented in Sections 3 and 4 of the reference paper. The
contents of the repository include:


Folder: Case_Cubic:                Corresponds to the cases with idealized homogeneous statistics (Sections 3).
             { main_GSI.m:         To construct the eigenmodes with GSI, GSI-P, and SFG for CFD computation of 
                                   homogeneous  turbulence (inhomogeneous boundaries can be present).
               main_TimeSeriesI.m: To generate temporal series at the inlet boundary of homogeneous turbulence using 
                                   eigenmodes constructed by GSI.
               f_TRL:              To construct the 3D covariance matrix of homogeneous turbulence with the TRL model. }
               Call orders of the main programs in Case_Cubic: main_GSI.m --> main_TimeSeries.m


Folder: Case_TBL:                  Corresponds to the turbulent boundary layer cases (Sections 4).
             { main_GSI.m:         To construct the eigenmodes with GSI, GSI-P, and SFG for CFD computation of TBL.
               main_TimeSeriesI.m: Generate temporal series at the inlet boundary of TBL using eigenmodes constructed by
                                   GSI.
               main_TRL.m:         To construct the 3D covariance matrix of TBL with the TRL model.}
               Call orders of the main programs in Case_TBTL: main_GSI.m --> main_TimeSeries.m


Folder: Lib_GSI:                   Contains the functions needed by the Cubic and TBL cases.
             { f_GSI_1D.m:         To optimize the phases of eigenmodes with the BFGS algorithm such that the integrated 
                                   variance of velocity divergence is minimized.
               f_J.m:              To calculate the integrated variance of velocity divergence as the cost function of the
                                   optimization problem.
               f_dJ.m:             To calculate the first derivative of the integrated variance of velocity divergence to the phases.
               f_Spectra_Pope.m:   Model spectrum presented in Turbulent Flows (Stephen B. Pope, 2000).
               f_Diff.m:           To calculate the differencing matrix for the first derivative. }
