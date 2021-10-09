function optimize_n_IR(dofdata,indata,LAMBDA)

% static=indata.static;
num_modes=indata.num_modes;
LAMBDA=LAMBDA(1:num_modes);

%mats=matfile('mats.mat','Writable',true);
mats=load('mats.mat');

if ~isfile('modes_I_global_store.mat')
    mats_I=struct;
    save('modes_I_global_store.mat','-struct','mats_I','-v7.3');
    mats_I=matfile('modes_I_global_store.mat','Writable',true);
    mats_I=store_modes_global(mats_I,mats,dofdata,indata);
    mats_I.Properties.Writable = false; % to prevent further changes
else
    mats_I=matfile('modes_I_global_store.mat','Writable',false);
end

mats_I=load('modes_I_global_store.mat');

indata.eigf.interface.method=1; % to keep modes based on a cutoff frequency

input_cma_opts
% Run CMAES to find the most probable values
[theta1,fval]=cmaes_v3_61('obj_n_IR',theta0,[],opts,mats,mats_I,dofdata,indata,LAMBDA);

%l=0;

% for k=1:1:20 % dummy loop to test speed
% 
%     indata.eigf.interface.target=k*4.6;
%     
%     [matdata,indata]=optmatassem_global([],mats,mats_I,dofdata,indata);
%        
%     if static==0
%         K_reduced=matdata.K_D_reduced;
%         M_reduced=matdata.M_D_reduced;
%     else
%         K_reduced=matdata.K_R_reduced;
%         M_reduced=matdata.M_R_reduced;
%     end
%     
%     LAMBDA_reduced=eigs(K_reduced,M_reduced,num_modes,'smallestabs');
%         
%     error=abs((sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)');
%     
%     l=l+1;
%     fprintf(['\nIteration #',num2str(l),'\n']);
%     %disp(error)
%     disp(indata.n_IR);
%     disp(max(error,[],'all'))
%     
% end



end