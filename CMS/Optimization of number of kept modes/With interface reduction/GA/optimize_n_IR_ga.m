clear

indata=load('indata.mat');
dofdata=load('dofdata.mat');
load LAMBDA.mat

% keep all modes up to "num_modes" to be used in error term
num_modes=indata.num_modes;
LAMBDA=LAMBDA(1:num_modes);

% file to keep progress of optimization
testfile=matfile('bestx_n_IR.mat','Writable',true);
testfile.bestx=[];

%mats=matfile('mats.mat','Writable',true);
mats=load('mats.mat');

% "mats_I" contains the stored interface modes
if ~isfile('modes_I_global_store.mat')
    mats_I=struct;
    save('modes_I_global_store.mat','-struct','mats_I','-v7.3');
    mats_I=matfile('modes_I_global_store.mat','Writable',true);
    mats_I=store_modes_global(mats_I,mats,dofdata,indata);
    mats_I.Properties.Writable = false; % to prevent further changes
    
    mats_I=load('modes_I_global_store.mat'); % load in RAM (faster)
else
    %mats_I=matfile('modes_I_global_store.mat','Writable',false);
    
    mats_I=load('modes_I_global_store.mat'); % load in RAM (faster)
end

% set to global interface reduction and no parametrization
indata.reduction_I=2;

% set number of kept modes explicitly using n_IR.
% This is necessary because here we are controlling n_IR directly
% and not the cuttoff frequency (as with CMA-ES).
indata.eigf.interface.method=0;

% the maximum allowable error (10^-2=1%)
tol=10^-2;

% number of design parameters -> number of interfaces=1
num_param=1;

% lower bounds: 1 mode per interface
lb=ones(num_param,1);

% upper bounds: the number of stored interface modes (in order to avoid
% solving the eigenproblem during optimization)
ub=length(mats_I.OMEGA_I_store);

% integer constraints: all design parameters are integers
intcons=1:num_param;

% this function writes the best individual in the file "testfile" at each
% iteration
param_test=@(options,state,flag) test_fun(options,state,flag,testfile);

% ga options
options=optimoptions('ga','Display','iter','plotfcn','gaplotbestindiv','OutputFcn',param_test);

% constraint function: must keep maximum error bellow "tol"
param_fun=@(theta) cons(theta,mats,mats_I,LAMBDA,tol,dofdata,indata);

% run the ga
tic
[x,fval,exitflag,output]=ga(@obj_fun,num_param,[],[],[],[],lb,ub,param_fun,intcons,options);
toc

testfile.Properties.Writable=false;
bestx=testfile.bestx;
figure
plot(bestx);


function f=obj_fun(theta)
    f=sum(theta.^2);
    %disp(theta)
end

function [c,ceq]=cons(theta,mats,mats_I,LAMBDA,tol,dofdata,indata) 
     
    indata.n_IR=theta;
    
    %tic
    [matdata,~,~]=optmatassem_global([],mats,mats_I,dofdata,indata);
    %toc
    
    %tic
    LAMBDA_reduced=eigs(matdata.K_D_reduced,matdata.M_D_reduced,length(LAMBDA),'smallestabs');
    %toc    
    
    error=max(abs((sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)'))./sqrt(LAMBDA)');
    
    c=(error-tol);
    fprintf('%50s --- %20s\n',num2str(theta),num2str(c));
    ceq=[];
end

function [state,options,optchanged] = test_fun(options,state,flag,testfile)
if flag=="iter"
    ibest = state.Best(end);
    ibest = find(state.Score == ibest,1,'last');
    bestx = state.Population(ibest,:);
    testfile.bestx=vertcat(testfile.bestx,bestx);
end
    optchanged=false;

end