%% ------------------------------------------------------------------------------------------ %%
%%       AAA      MM       MM      AAA      LL        GGGGGGGGGGG      AAA      MM       MM   %%
%%      AA AA     MMM     MMM     AA AA     LL        GGG     GGG     AA AA     MMM     MMM   %%
%%      AA AA     MMMM   MMMM     AA AA     LL        GG       GG     AA AA     MMMM   MMMM   %%
%%     AA   AA    MM MM MM MM    AA   AA    LL        GGG     GGG    AA   AA    MM MM MM MM   %%
%%    AAAAAAAAA   MM  MMM  MM   AAAAAAAAA   LL        GGGGGGGGGGG   AAAAAAAAA   MM  MMM  MM   %%
%%    AA     AA   MM       MM   AA     AA   LL                 GG   AA     AA   MM       MM   %%
%%   AA       AA  MM       MM  AA       AA  LLLLLLLL           GG  AA       AA  MM       MM   %%
%%   AA       AA  MM       MM  AA       AA  LLLLLLLL  GGGGGGGGGGG  AA       AA  MM       MM   %%
%% ------------------------------------------------------------------------------------------ %%
%% ------------------ The AMALGAM multiobjective optimization algorithm --------------------- %%
%%                                                                                            %%
%% This general purpose MATLAB code is designed to find a set of parameter values that        %%
%% defines the Pareto trade-off surface corresponding to a vector of different objective      %%  
%% functions. In principle, each Pareto solution is a different weighting of the objectives   %% 
%% used. Therefore, one could use multiple trials with a single objective optimization        %% 
%% algorithms using diferent values of the weights to find different Pareto solutions.        %% 
%% However, various contributions to the optimization literature have demonstrated that this  %% 
%% approach is rather inefficient. The AMALGAM code developed herein is designed to find an   %% 
%% approximation of the Pareto solution set within a single optimization run. The AMALGAM     %% 
%% method combines two new concepts, simultaneous multimethod search, and self-adaptive       %% 
%% offspring creation, to ensure a fast, reliable, and computationally efficient solution to  %% 
%% multiobjective optimization problems. This method is called a multi-algorithm, genetically %% 
%% adaptive multiobjective, or AMALGAM, method, to evoke the image of a procedure that blends %% 
%% the attributes of the best available individual optimization algorithms.                   %%                                         
%%                                                                                            %%
%% ------------------------------------------------------------------------------------------ %%
%%                                                                                            %%
%% SYNOPSIS:                                                                                  %%
%%                                                                                            %%
%%        [X,F,output,Z,sim] = AMALGAM(AMALGAMPar,Func_name,Par_info);                        %%
%%        [X,F,output,Z,sim] = AMALGAM(AMALGAMPar,Func_name,Par_info,options);                %%
%%        [X,F,output,Z,sim] = AMALGAM(AMALGAMPar,Func_name,Par_info,options,func_in);        %%
%%        [X,F,output,Z,sim] = AMALGAM(AMALGAMPar,Func_name,Par_info,options,func_in,Fpar);   %%
%%                                                                                            %%
%% Input:    AMALGAMPar = structure with AMALGAM settings/parameters                          %%
%%           Func_name = name of the function or model that returns objective functions       %%
%%           Par_info = optional structure with parameter ranges                              %%
%%           Fpareto = optional vector with Pareto solution set (benchmark problems)          %%
%%           options = optional structure with additional settings                            %%
%%           func_in = optional variable that user can pass to function                       %%
%%                                                                                            %%
%% Output:   X = final population (matrix)                                                    %%
%%           F = final objective function values of "X" (matrix)                              %%
%%           output = structure with several output arguments computed by AMALGAM (structure) %%
%%           Z = archive of all past populations augmented with X (matrix)                    %%
%%           sim (optional) = Model simulations of Pareto solutions (see example 6 and 7)     %%
%%                                                                                            %%
%% ------------------------------------------------------------------------------------------ %%
%%                                                                                            %%
%% This algorithm has been described in:                                                      %%
%%                                                                                            %%
%% Vrugt, J.A., and B.A. Robinson, Improved evolutionary optimization from genetically        %%
%%    adaptive multimethod search, Proceedings of the National Academy of Sciences of the     %%
%%    United States of America, 104, 708 - 711, doi:10.1073/pnas.0610471104, 2007.            %%
%% Vrugt, J.A., B.A. Robinson, and J.M. Hyman, Self-adaptive multimethod search for global    %%
%%    optimization in real-parameter spaces, IEEE Transactions on Evolutionary Computation,   %%
%%    13(2), 243-259, doi:10.1109/TEVC.2008.924428, 2009.			                          %%
%%                                                                                            %%
%% For more information please read:                                                          %%
%%                                                                                            %%
%% Vrugt, J.A., H.V. Gupta, L.A. Bastidas, W. Bouten, and S. Sorooshian, Effective and        %%
%%    efficient algorithm for multi-objective optimization of hydrologic models, Water        %% 
%%    Resources Research, 39(8), art. No. 1214, doi:10.1029/2002WR001746, 2003.               %%                        
%% Schoups, G.H., J.W. Hopmans, C.A. Young, J.A. Vrugt, and W.W.Wallender, Multi-objective    %%
%%    optimization of a regional spatially-distributed subsurface water flow model, Journal   %%
%%    of Hydrology, 20 - 48, 311(1-4), doi:10.1016/j.jhydrol.2005.01.001, 2005.               %%
%% Vrugt, J.A., P.H. Stauffer, T. Wöhling, B.A. Robinson, and V.V. Vesselinov, Inverse        %% 
%%    modeling of subsurface flow and transport properties: A review with new developments,   %%
%%    Vadose Zone Journal, 7(2), 843 - 864, doi:10.2136/vzj2007.0078, 2008.                   %%
%% Wöhling, T., J.A. Vrugt, and G.F. Barkle, Comparison of three multiobjective optimization  %%
%%    algorithms for inverse modeling of vadose zone hydraulic properties, Soil Science       %% 
%%    Society of America Journal, 72, 305 - 319, doi:10.2136/sssaj2007.0176, 2008.            %%
%% Wöhling, T., and J.A. Vrugt, Combining multi-objective optimization and Bayesian model     %%
%%    averaging to calibrate forecast ensembles of soil hydraulic models, Water Resources     %%
%%    Research, 44, W12432, doi:10.1029/2008WR007154, 2008.                                   %%
%%                                                                                            %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                                            %%
%% AMALGAM code developed by Jasper A. Vrugt, University of California Irvine: jasper@uci.edu %%
%%                                                                                            %%
%% Version 0.5:   June 2006                                                                   %%
%% Version 1.0:   January 2009    Cleaned up source code and implemented 4 test problems      %%
%% Version 1.1:   January 2010    Flexible population size and no need divide by # algorithms %%
%% Version 1.2:   August 2010     Sampling from prior distribution                            %%
%% Version 1.3:   May 2014        Varous updates - cleaning and improved speed ranking        %%
%% Version 1.4:   Januari 2014    Parallellization using parfor (done if CPU > 1)             %%
%% Version 2.0:   May 2017        Final implementation (many updates to code)                 %%
%%                December 2017   Discrete sampling ( aka DREAM_D )		              %%		
%%                                                                                            %%
%% ------------------------------------------------------------------------------------------ %%

%% ------------------------------------------------------------------------------------------ %%
%%                                                                                            %%
%%                       PLEASE CHECK MANUAL OF AMALGAM (ON MY WEBSITE)                       %% 
%%                                                                                            %%
%% Vrugt, J.A., Multi-criteria optimization using the AMALGAM software package: Theory,       %%
%%    concepts, and MATLAB implementation, Manual, p. 1-70, 2015.                             %%
%%                                                                                            %%
%% ------------------------------------------------------------------------------------------ %%

%% ------------------------------------------------------------------------------------------ %%
%%                                                                                            %%
%%   NOTE: EXPLICIT PRIOR SAMPLING DISTRIBUTIONS CAN BE USED IN AMALGAM: CHECK DREAM MANUAL   %%
%%                                                                                            %%
%% Vrugt, J.A., Markov chain Monte Carlo simulation using the DREAM software package: Theory, %%
%%    concepts, and MATLAB Implementation, Environmental Modelling & Software, 75, 273-316,   %%
%%    doi:10.1016/j.envsoft.2015.08.013.                                                      %%
%%                                                                                            %%
%% ------------------------------------------------------------------------------------------ %%

%% Check:  http://faculty.sites.uci.edu/jasper
%% Papers: http://faculty.sites.uci.edu/jasper/publications/
%% Google Scholar: https://scholar.google.com/citations?user=zkNXecUAAAAJ&hl=nl

%% Different test examples (1 - 5 are synthetic problems, 6 is a real world problem)
%% example 1: ZDT1: test function
%% example 2: ZDT2: test function
%% example 3: ZDT3: test function
%% example 4: ZDT4: test function
%% example 5: ZDT6: test function
%% example 6: real-world example using streamflow simulation with hmodel
%% example 7: watershed modeling with HYMOD using driven and nondriven part hydrograph
%% example 8: Bayesian model averaging using RMSE, IS and CRPS as metrics
%% example 9: Multicriteria BMA model calibration using temperature ensemble 
%% example 10: Multicriteria BMA model calibration using sea-level pressure ensemble 
%% example 11: ZDT4: test function but with discrete (e.g. integer) sampling

%% ########################################################################
%%   Func_name: Name of the function script of the model/function
%% ########################################################################
%%                        CASE STUDY DEPENDENT
%% ------------------------------------------------------------------------
%% Func_name                 % Name of the model function script (.m file)
%% ------------------------------------------------------------------------

%% ########################################################################
%%   AMALGAMPar: Computational setup AMALGAM and algorithmic parameters
%% ########################################################################
%%                         CASE STUDY DEPENDENT
%% ------------------------------------------------------------------------
%% AMALGAMPar.d              % Dimensionality Pareto distribution
%% AMALGAMPar.N              % Population size
%% AMALGAMPar.T              % Number of generations?
%% AMALGAMPar.m              % Number of objective functions?
%% ------------------------------------------------------------------------
%%                           DEFAULT VALUES
%% ------------------------------------------------------------------------
%% AMALGAMPar.rec_methods   % Recombination methods  : {'GA','PSO','AMS','DE'}
%% AMALGAMPar.beta_1        % DE scaling factor      : @(N) unifrnd(0.6,1.0,N,1)
%% AMALGAMPar.beta_2        % DE scaling factor      : @(N) unifrnd(0.2,0.6,N,1)
%% AMALGAMPar.c_1           % PSO social factor      : 1.5
%% AMALGAMPar.c_2           % PSO cognitive factor   : 1.5
%% AMALGAMPar.varphi        % PSO inertia factor     : @(N) unifrnd(0.5,1.0,N,1)
%% AMALGAMPar.p_CR          % NSGA-II crossover rate : 0.9
%% AMALGAMPar.p_M           % NSGA_II mutation rate  : 1/d
%% AMALGAMPar.eta_C         % NSGA-II mutation index : 10
%% AMALGAMPar.eta_M         % NSGA-II mutation index : 50
%% AMALGAMPar.gamma         % AMS jump rate          : (2.38/sqrt(d))^2
%% AMALGAMPar.K             % Thinning rate          : 1 (no thinning of Z)
%% AMALGAMPar.p_min         % Min. selection prob.   : 0.05
%% ------------------------------------------------------------------------

%% ########################################################################
%%   Par_info: All information about the parameter space and prior
%% ########################################################################
%%                        CASE STUDY DEPENDENT
%% ------------------------------------------------------------------------
%% Par_info.initial          % Initial sampling distribution ('uniform'/'latin'/'normal'/'prior')
%% Par_info.min              % If 'latin', min parameter values
%% Par_info.max              % If 'latin', max parameter values
%% Par_info.prior            % Marginal prior distribution of each parameter
%% Par_info.mu               % If 'normal', mean of initial parameter values
%% Par_info.cov              % If 'normal', covariance matrix parameters
%% Par_info.boundhandling    % Boundary handling ('reflect','bound','fold')
%% Par_info.steps            % Number of intervals each parameter (discrete sampling)
%% ------------------------------------------------------------------------
%%                          DEFAULT VALUES
%% ------------------------------------------------------------------------
%% Par_info.boundhandling = 'none'   % no boundary handling (unbounded problem)
%% ------------------------------------------------------------------------

%% ########################################################################
%%   options: Structure with optional settings
%% ########################################################################
%%                             OPTIONAL
%% ------------------------------------------------------------------------
%% options.parallel          % Multi-core computation chains?
%% options.IO                % If parallel, IO writing model?
%% options.save              % Save DREAM output during the run?
%% options.ranking           % Pareto Ranking code, 'MATLAB' (default) or 'C' (faster)
%% options.density           % Which density of points 'crowding' (default) or 'strength'
%% options.modout            % Return model simulatons? 'no' (default) or 'yes'
%% options.restart           % Restart run ( continue previous run - options.save must be 'yes')
%% options.print             % Print to screen tables/figures (postprocessor)
%% ------------------------------------------------------------------------

%% ########################################################################
%%  Fpareto: Matrix ( Npar x d ) with Pareto solutions (synthetic problems)
%%  NOTE: Existing IGD.mexw64 in zip file compiled for 64 bit machine
%%  NOTE: If this gives error recompile IGD.cpp ("mex IGD.cpp")
%%  NOTE: If you do not have mex compiler and IGD gives errors just specify
%%  NOTE: Fpar = []; 
%% ########################################################################
%%                             OPTIONAL
%% -------------------------------------------------------------------------
%% Fpar = []                % Empty -> no Pareto solutions specified
%% -------------------------------------------------------------------------

% Add main AMALGAM directory and underlying postprocessing directory
addpath(pwd,[pwd,'/postprocessing']);
% Now go to example 1
cd example_1
% And you can now execute example_1 using (uncomment)
% example_1