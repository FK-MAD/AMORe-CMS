%% TO DO
%1) Create all associated matrices for model without soil with COMSOL at home
%2) Copy the folders on the big computer
%3) Choose r (by hand)
%4) Run optimization for n_IR (one time) and n_IR_l (multiple times)->Choose n_IR and n_IR_l
%5) Create parametrized matrices in big computer
%6) Choose 10-20 random sample points and calculate LAMBDA of unreduced (from COMSOL) there
%7) For the same sample points, check the accuracy of each parametrization (play with scatter)
%8) Give parametrized matrices to Costas to run them in Bayesian framework

%% input

%group_l={[1:3],[4:7],8,[9:11],12,13,14,[15:18],19,[20:22]}

% ---parallelization settings---
indata.num_workers=6; % number of workers to use
% ---


% ---optimization settings---
indata.num_modes=20; % number of first modes to use in error term between unreduced and reduced model
% ---


% ---matrix assembly settings---
% method of building the block-diagonal matrices. It might help with very
% big matrices that may not fit in RAM.
% 1 -> use files mats_S_k.mat, read variables in a for-loop and build the block-diagonal incrementally
% 2 -> use files mats_S_k.mat, read variables in a for-loop and build the block-diagonal once
% 3 -> use file mats_S.mat and build the block-diagonal once
indata.blkdiag_method=2;
% ---


% ---reduction method---
% without parametrization: 1=no interface reduction | 2=global interface reduction | 3=local interface reduction
% with parametrization:    4=no interface reduction | 5=global interface reduction | 6=local interface reduction
indata.reduction_I=1;
% ---


% ---use of static correction---
% 0=without static correction | 1=with static correction
indata.static=0;
% ---


% ---kept modes for component groups---
% method of calculating the kept modes
% 0=explicitly using n_id_S | 1=until the target eigenfrequency for each group is reached
indata.eigf.group.method=1; 

% all vectors have:
% rows=1
% columns>=number of component groups (will run normally if more columns than component groups exist)

% this is used if method=0
indata.n_id_S=50*ones(1,100); % kept fixed-interface normal modes for each group of components

% this is used if method=1
r=8*ones(1,100); % multiplication constant used to define the target frequency
indata.eigf.group.multiplier=r;
indata.eigf.group.target=r*4.5; % target eigenfrequency (Hz) for each group of components

% this controls the way modes are searched
indata.eigf.group.max=500*ones(1,100); % maximum allowed number of modes
indata.eigf.group.step=50*ones(1,100); % increase in the number of calculated modes if target is not reached
indata.eigf.group.init=50*ones(1,100); % initial number of calculated modes
% ---


% ---stored modes for component groups---
% they are computed once and used when updating matrices during
% optimization of r. They should be enough to avoid solving the eigenproblem
% during optimization.
% If you don't want any stored modes select:
% indata.eigf.group.method_store=0 and
% indata.n_id_S_store=0*ones(1,100)

% everything here works similarly to kept modes (same logic)

indata.eigf.group.method_store=1; % 0 or 1

% this is used if method_store=0
indata.n_id_S_store=50*ones(1,100); % large values

% this is used if method_store=1
indata.eigf.group.target_store=20*4.5*ones(1,100); % large cutoff frequency

% this controls the way modes are searched
indata.eigf.group.max_store=500*ones(1,100); % large values
indata.eigf.group.step_store=50*ones(1,100); % large values
indata.eigf.group.init_store=50*ones(1,100); % large values
% ---


% ---kept modes for interfaces---
% method of calculating the kept modes
% 0=explicitly using n_IR (for global reduction) or n_IR_l (for local reduction) | 1=until the target eigenfrequency for each interface is reached
indata.eigf.interface.method=1;

% all vectors have:
% rows=1
% columns>=number of interfaces (will run normally if more columns than interfaces exist)

% if global interface reduction is selected:
% only the first element of the target, max, step and init vectors is used
% (there is only one interface)

% this is used if method=0 and global reduction is selected
indata.n_IR=50; % kept interface modes for all interfaces (global reduction)

% this is used if method=0 and local reduction is selected
indata.n_IR_l=50*ones(1,100); % kept interface modes for each interface (local reduction)

% this is used if method=1
v=3*ones(1,100); % multiplication constant used to define the target frequency
% ~1% error-> metsovo: 3 for global=70 for local
indata.eigf.interface.multiplier=v;
indata.eigf.interface.target=v*4.5; % target eigenfrequency (Hz) for each interface

% this controls the way modes are searched
indata.eigf.interface.max=500*ones(1,100); % maximum allowed number of modes
indata.eigf.interface.step=50*ones(1,100); % increase in the number of calculated modes if target is not reached
indata.eigf.interface.init=50*ones(1,100); % initial number of calculated modes
% ---


% ---stored modes for interfaces---
% they are computed once and used when updating matrices during
% optimization of v. They should be enough to avoid solving the eigenproblem
% during optimization.
% If you don't want any stored modes select:
% indata.eigf.interface.method_store=0 and indata.n_IR_l_store=0*ones(1,100)

% everything here works similarly to kept modes (same logic)

indata.eigf.interface.method_store=1; % 0 or 1

% this is used if method_store=0 and global reduction is selected
indata.n_IR_store=50; % for global reduction, large value

% this is used if method_store=0 and local reduction is selected
indata.n_IR_l_store=50*ones(1,100); % for local reduction, large values

% this is used if method_store=1
indata.eigf.interface.target_store=70*4.5*ones(1,100); % large cutoff frequency

% this controls the way modes are searched
indata.eigf.interface.max_store=500*ones(1,100); % large values
indata.eigf.interface.step_store=50*ones(1,100); % large values
indata.eigf.interface.init_store=50*ones(1,100); % large values
% ---


% ---material properties---
% all vectors have:
% rows>=number of component groups (will run normally if more columns than component groups exist)
% columns=1

indata.E=37*10^9*ones(17,1); % Young's modulus [Pa] for each group of components. most groups are deck components -> 37 GPa
indata.E([6,11,14])=34*10^9; % groups 6, 11 and 14 are piers -> 34 GPa

indata.nu=.2*ones(17,1); % Poisson's ratio

indata.rho=2548*ones(17,1); % density [kg/m^3]
% ---


% ---general parametrization settings---
% interpolation scheme used in interpolation of interface modes (global or local reduction)
indata.quad_interp=0; % 0=linear interpolation | 1=quadratic interpolation

% static correction method
indata.invariant=0; % 0=full static correction | 1=invariant assumption


% functions of the model parameters -> one entry for each model parameter
% func_g applies on mass matrix (see Eq. (2.3))
indata.func_g=repmat({@(x) 1},1,15);

% func_h applies on stiffnes matrix (see Eq. (2.4))
indata.func_h=repmat({@(x) x},1,15);


% vectors have:
% rows=number of model parameters
% columns=1

% sample point where reduced matrices are calculated [1.1 0.9 0.85 1.15]'.
% This is used to test the code. Normally, every sample point is generated
% during the stochastic simulation process.
indata.theta_k=[1 1 1 1 1 1 1 1 1 1]';

% nominal point used in the invariant assumption of static correction. See
% page 42.
indata.theta_nom=ones(15,1);

% nominal point used in interpolation of interface modes (global or local
% reduction). See page 50.
indata.theta_0=ones(15,1);
% ---


% ---settings concerning support points (if parametrization is used)---
% used in interpolation of interface modes (global or local reduction)

% 'scatter_theta_l' is the fraction of theta_0 that the support points are
% scattered around theta_0 (can be different for each parameter)
% e.g. for a model with 2 parameters and theta_0=[1;1]:
% scatter_theta_l(1)=.1 -> support of parameter 1=[.9,1.1]
% scatter_theta_l(2)=.2 -> support of parameter 2=[.8,1.2]
indata.scatter_theta_l=.6*ones(length(indata.theta_0),1);

% 'simplex' -> for dimension n there are needed n+1 vertices to create the convex hull
% 1D -> 2 points -> line segment
% 2D -> 3 points -> triangle
% 3D -> 4 points -> tetrahedron
% 4D -> 5 points -> 5-cell (4-simplex)
% ...

% 'hypercube' -> for dimension n there are needed 2^n vertices to create the convex hull
% 1D -> 2 points -> line segment
% 2D -> 4 points -> square
% 3D -> 8 points -> cube
% 4D -> 16 points -> 4-cube (hypercube,tesseract)
% ...

% both methods provide as few support points as possible using smart 
% merging rules
indata.method_theta_l='simplex'; % 'simplex' or 'hypercube'

%indata.theta_l=mvnrnd(indata.theta_0,.02*eye(length(indata.theta_0)),10)'; % random sampling
%indata.theta_l=lhsnorm(indata.theta_k,.1*eye(length(indata.theta_0)),20)'; % latin hypercube sampling
%theta_l=[1.2 0.8 0.8 0.8;0.8 1.2 0.8 0.8;0.8 0.8 1.2 0.8;0.8 0.8 0.8 1.2]';
% ---


%% pass additional input data to structure "indata"

indata.filename=filename;
indata.save_dir=save_dir;

indata.S_0=S_0; % groups that are independent of model parameters
indata.S_j=S_j; % groups that depend of model parameters. cell 1,2,... contains the groups that depend on parameter 1,2,...
indata.n_theta=length(indata.S_j); % number of parameters. n_theta=length(func_g);
%indata.L=size(indata.theta_l,2); % number of support points

indata.group_S=group_S; % grouping of geometrical domains. cell 1,2,... contains the domains that make up group 1,2,...
indata.n_id=sum(indata.n_id_S); % total number of kept fixed-interface normal modes
indata.N_S=length(indata.group_S); % number of groups of components
indata.n_IRL=sum(indata.n_IR_l); % total number of kept interface modes using local reduction

indata.n_DIL=indata.n_id+indata.n_IRL; % dimension of reduced matrices if local reduction is used
indata.n_DI=indata.n_id+indata.n_IR; % dimension of reduced matrices if global reduction is used