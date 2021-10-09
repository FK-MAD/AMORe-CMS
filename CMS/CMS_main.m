%close all
clear

% ---model name---
filename='Metsovo_normal_mesh_22para'; % the filename of the model to load
% ---

% ---save directory---
save_dir='D:\Philip\Desktop\Metsovo - normal mesh - 22 parameters'; % folder to save all files
% ---

%% Pre-processing

% Add path of COMSOL models and needed functions to the current folder
addpath(genpath('../COMSOL/'));
addpath(genpath('../CMS/'));

% the import statements below make all model and model utility methods available
import com.comsol.model.*;
import com.comsol.model.util.*;

%ModelUtil.clear(); % clear all models loaded in COMSOL server
%mphlaunch(filename); % launces the model in COMSOL

model_tags=mphtags; % tags of models already loaded in COMSOL server
if ismember(filename,string(model_tags(:))) % if model is already loaded in COMSOL server
    model=ModelUtil.model(filename); % import it from server
else
    model=mphload(filename,filename); % load it and create a linked model object in MATLAB
end

[group_S,S_0,S_j]=groupassem(model); % extract the dependency of component groups on model parameters from COMSOL

CMS_input_metsovo_22para; % run input file

current=cd(save_dir); % go to specified folder


%% Construct needed matrices

% this checks for moddata.mat -> contains relevant information about the
% model. Needs COMSOL server.
if ~isfile('moddata.mat') % file does not exist
    moddata=modeldata(model,indata);
    save('moddata.mat','-struct','moddata','-v7.3');
else % file exists
    fprintf('\nmoddata.mat already exists and will be used instead of creating a new one...\n\n');
    moddata=load('moddata.mat');
end

% this checks for dofdata.mat -> contains information about the DOFs of the
% model. Needs COMSOL server.
if ~isfile('dofdata.mat') % file does not exist
    [dofdata,indata]=dofassem(model,moddata,indata);
    save('dofdata.mat','-struct','dofdata','-v7.3');
else % file exists
    fprintf('dofdata.mat already exists and will be used instead of creating a new one...\n\n');
    dofdata=load('dofdata.mat');
    
    % use the interface selection ("group_l") contained in "indata.mat"
    % which corresponds to "dofdata.mat" that already exists 
    load('indata.mat','group_l');
    indata.group_l=group_l;
    clear group_l
end

% moddata might not be needed. Consider clearing the variable to save RAM.
% clear moddata

% this generates support points based on input. It needs to run after the creation
% of "dofdata.mat"  because it uses information contained there
tic
indata=make_support_points(indata,dofdata); 
elapsed_time_support=toc;

% this makes files "mats_S_k.mat" (k=1,2,...,N_S) where the partitioned mass and stiffness matrices
% of each component group are stored. More matrices will be added in these
% files later. This is the last function that needs the COMSOL server.
if ~isfile('mats_S_1.mat') % if this file doesn't exist, I assume that the corresponding files for all other groups don't exist either
    matextract_COMSOL(model,dofdata,indata)
end  

% this checks input concerning kept component modes (interface modes). If
% vectors contain more elements than component groups (interfaces) there will
% be a warning. If vectors contain less elements there will be an error.
indata=error_check(dofdata,indata);

% save indata.mat
save('indata.mat','-struct','indata','-v7.3');

% visualize geometry and DOFs
%correspondence_I=geomviz(model,dofdata,indata,'','');
%dofviz(moddata,dofdata,indata,'all');

% this creates a parallel pool with selected input. If running pool has
% different number of workers than input, the pool is shut down and a new
% pool is started.
pool=gcp('nocreate');
if isempty(pool) ||  pool.NumWorkers~=indata.num_workers
    delete(pool);
    parpool('local',indata.num_workers); % use specified number of cores
end

% this finalizes "mats_S_k.mat" (k=1,2,...,N_S) and creates "mats_S.mat"
% and "mats.mat". If they exist, it performs checks to see if they
% correspond to user input. It can update them to match the input, if the
% user wishes it.
indata=submatassem(dofdata,indata);

% this creates the final reduced matrices. Edit the function "matassem.m"
% to control if matrices are stored in RAM (faster) or in Matfiles (slower
% but needs less RAM).
[matdata,indata]=matassem(dofdata,indata);

% model updating. These functions are used during a stochastic simulation
% process. They should be modified accordingly.
% if indata.reduction_I==4
%     [matdata,indata]=model_updating_unreduced(dofdata,indata);
% elseif indata.reduction_I==5
%     [matdata,indata]=model_updating_global(dofdata,indata);
% elseif indata.reduction_I==6
%     [matdata,indata]=model_updating_local(dofdata,indata);
% end

%% Visualize error between eigenvalues and modeshapes (MAC) 

% visualize error between eigenvalues of reduced and unreduced model
LAMBDA=errorviz(indata,matdata);

% MAC visualization. Use this if transformation matrix "T_D" is included
% in "matdata". Normally this transformation matrix is not included because
% it requires additional calculations.
%
% T_G=dofdata.T_G;
% T_D=matdata.T_D_reduced;
% 
% V_D_transformed=T_G*T_D*V_D;
% MACviz(V,V_D_transformed,modes)

cd(current);