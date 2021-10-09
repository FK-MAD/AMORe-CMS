clear

% load all necessary variables in RAM (for speed)
mats=load('mats.mat');
indata=load('indata.mat');
dofdata=load('dofdata.mat');
load('mats_S.mat','PSI_ib_S');
load('mats_S.mat','M_ii_S');
load('mats_S.mat','M_ib_S');
load('mats_S.mat','K_ii_S');
load('mats_S.mat','PHI_id_S_store');
load('mats_S.mat','LAMBDA_id_S_store');
load LAMBDA.mat

% keep all modes up to "num_modes" to be used in error term
num_modes=indata.num_modes;
LAMBDA=LAMBDA(1:num_modes);

% file to keep progress of optimization
testfile=matfile('bestx_n_id_S.mat','Writable',true);
testfile.bestx=[];

% set to no interface reduction and no parametrization
indata.reduction_I=1;

% set number of kept modes explicitly using n_id_S.
% This is necessary because here we are controlling n_id_S directly
% and not the cuttoff frequency (as with CMA-ES).
indata.eigf.group.method=0;

% the maximum allowable error (10^-2=1%)
tol=10^-2;

% number of design parameters -> number of component groups
num_param=indata.N_S;

% lower bounds: 1 mode per component group
lb=ones(num_param,1);

% upper bounds: the number of stored modes for each component group (in order to avoid
% solving the eigenproblem during optimization)
[ub,~]=cellfun(@size,LAMBDA_id_S_store');

% integer constraints: all design parameters are integers
intcons=1:num_param;

% this function writes the best individual in the file "testfile" at each
% iteration
param_test=@(options,state,flag) test_fun(options,state,flag,testfile);

% ga options
options=optimoptions('ga','Display','iter','plotfcn','gaplotbestindiv','OutputFcn',param_test);

% constraint function: must keep maximum error bellow "tol"
param_fun=@(theta) cons(theta,mats,PSI_ib_S,M_ii_S,M_ib_S,K_ii_S,...
    PHI_id_S_store,LAMBDA_id_S_store,LAMBDA,tol,dofdata,indata);

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

function [c,ceq]=cons(theta,mats,PSI_ib_S,M_ii_S,M_ib_S,K_ii_S,...
    PHI_id_S_store,LAMBDA_id_S_store,LAMBDA,tol,dofdata,indata)

    indata.n_id_S=theta;
    
    %tic
    [mats,~,~,indata]=...
        update_mats_RAM(mats,PSI_ib_S,M_ii_S,M_ib_S,...
        K_ii_S,PHI_id_S_store,LAMBDA_id_S_store,dofdata,indata);
    %toc

    %tic
    matdata=matassem_unreduced([],mats,dofdata,indata);
    %toc

    %tic
    LAMBDA_reduced=eigs(matdata.K_D_reduced,matdata.M_D_reduced,length(LAMBDA),'smallestabs');
    %toc

    error=max(abs(sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)');

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