% INPUT_CMA_OPTS defines the input variables used in CMA-ES
% Set the options for CMA-ES
%--------------------------------------------------------------------------
% IMPORTANT CMA-ES parameters
%--------------------------------------------------------------------------
% Max number of generations
opts.MaxFunEvals=50000;
% PopSize = Number of samples per generation. Important parameter
% Recommended from Hansen as minimum : d * sqrt (d) d = parameter dimension
opts.PopSize=100; 
% Lower and upper bound of the design variables x
%opts.LBounds = [x_1_Lower  x_2_Lower]';
%opts.UBounds = [x_1_Upper  x_2_Upper]';
opts.LBounds =.01*ones(indata.N_S,1);
opts.UBounds =10*ones(indata.N_S,1);
%--------------------------------------------------------------------------
% Initial values of the design variables
% theta0 = [theta0(1) theta0(2)]';
theta0=7*ones(indata.N_S,1);
% Plot option
opts.LogPlot=1;
%--------------------------------------------------------------------------
% No need to change the following options
opts.CMA.active = 1;
opts.Resume=0;
opts.Noise.on = 0;
opts.EvalInitialx='no';
opts.EvalParallel='no';
opts.LogModulo=1;
%--------------------------------------------------------------------------
