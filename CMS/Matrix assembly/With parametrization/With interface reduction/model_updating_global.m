function [matdata,indata]=model_updating_global(dofdata,indata)

if ~isfile('submatdata_M_I_K_I.mat')
    submatdata_M_I_K_I=matfile('submatdata_M_I_K_I.mat','Writable',true);
    mats_S=matfile('mats_S.mat','Writable',false);
    submatdata_M_I_K_I=submatassem_M_I_K_I(submatdata_M_I_K_I,mats_S,dofdata,indata);
    submatdata_M_I_K_I.Properties.Writable = false; % to prevent further changes
end
submatdata_M_I_K_I=load('submatdata_M_I_K_I.mat');

if ~isfile('submatdata_interp_global.mat')
    submatdata_interp=matfile('submatdata_interp_global.mat','Writable',true);
    [submatdata_interp,indata]=submatassem_interp(submatdata_interp,submatdata_M_I_K_I,dofdata,indata);
    submatdata_interp.Properties.Writable = false; % to prevent further changes
    %save('indata.mat','-struct','indata','-v7.3');
end
submatdata_interp=load('submatdata_interp_global.mat');

n_IR_correct=size(submatdata_interp.YPSILON_I_theta_0,2);
if ~isequal(n_IR_correct,indata.n_IR)
    warning(['Input vector n_IR does not correspond to the number of kept interface modes.'...
        ' Will change the input vector to correspond properly.']);
    indata.n_IR=n_IR_correct;
    indata.n_DI=indata.n_id+indata.n_IR;
end

% % to ensure that the number of kept interface modes corresponds to the
% % input file (submatdata_interp). Else, there is error when calculating
% % interface modes in a sample point which does not belong in the convex
% % hull of the existing support points.
% indata.eigf.interface.method=0;

if ~isfile('submatdata_5.mat') 
    submatdata=matfile('submatdata_5.mat','Writable',true);
    
    mats=matfile('mats.mat','Writable',false);
    mats_S=matfile('mats_S.mat','Writable',false);
    
    indata.theta_k=indata.theta_nom;
    mats_interp_nom=matassem_interp([],submatdata_interp,submatdata_M_I_K_I,dofdata,indata);
    
    submatdata=submatassem_para_global(submatdata,mats_S,mats,mats_interp_nom,dofdata,indata);
    submatdata.Properties.Writable = false; % to prevent further changes
end
submatdata=load('submatdata_5.mat');

tic   
[mats_interp_k,submatdata_interp,indata]=interpolate_I([],submatdata_interp,submatdata_M_I_K_I,dofdata,indata);
toc

tic
matdata=matassem_para_global([],submatdata,mats_interp_k,dofdata,indata);
toc
    
% l=0;
% for k=1:10 % dummy loop to test speed
%     l=l+1;
%     fprintf(['\nIteration #',num2str(l),'\n']);
%     indata.theta_k=normrnd(1,.05,4,1); % change the sample point      
%     
%     tic   
%     [mats_interp_k,submatdata_interp,indata]=interpolate_I([],submatdata_interp,submatdata_M_I_K_I,dofdata,indata);
%     toc
%     
%     tic
%     matdata=matassem_para_global([],submatdata,mats_interp_k,dofdata,indata);
%     toc
% end

end