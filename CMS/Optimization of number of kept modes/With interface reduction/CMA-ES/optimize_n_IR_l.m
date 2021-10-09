function optimize_n_IR_l(dofdata,indata,LAMBDA)

% N_I=dofdata.N_I;
% static=indata.static;
num_modes=indata.num_modes;
LAMBDA=LAMBDA(1:num_modes);

%mats=matfile('mats.mat','Writable',true);
mats=load('mats.mat');

if ~isfile('modes_I_local_store.mat')
    mats_I=struct;
    save('modes_I_local_store.mat','-struct','mats_I','-v7.3');
    mats_I=matfile('modes_I_local_store.mat','Writable',true);
    mats_I=store_modes_local(mats_I,mats,dofdata,indata);
    mats_I.Properties.Writable = false; % to prevent further changes
else
    mats_I=matfile('modes_I_local_store.mat','Writable',false);
end

mats_I=load('modes_I_local_store.mat');

indata.eigf.interface.method=1; % to keep modes based on a cutoff frequency

input_cma_opts
% Run CMAES to find the most probable values
[theta1,fval]=cmaes_v3_61('obj_n_IR_l',theta0,[],opts,mats,mats_I,dofdata,indata,LAMBDA);


% l=0;
% 
% for k=10:10:100 % dummy loop to test speed
% 
%     indata.eigf.interface.target=k*4.6*ones(1,N_I);
%     
%     % to increase speed consider calculating K_II and M_II as full (if
%     % memory allows for this) in optmatassem_local
%     % K_II=full(YPSILON_Ill)'*full(K_I)*full(YPSILON_Ill);
%     % M_II=full(YPSILON_Ill)'*full(M_I)*full(YPSILON_Ill);
%     tic
%     [matdata,indata]=optmatassem_local([],mats,mats_I,dofdata,indata);
%     toc
%     
%     if static==0
%         K_reduced=matdata.K_D_reduced;
%         M_reduced=matdata.M_D_reduced;
%     else
%         K_reduced=matdata.K_R_reduced;
%         M_reduced=matdata.M_R_reduced;
%     end
%     
%     tic
%     LAMBDA_reduced=eigs(K_reduced,M_reduced,num_modes,'smallestabs');
%     toc    
%     
%     error=abs((sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)')*100;
%     
%     l=l+1;
%     fprintf(['\nIteration #',num2str(l),'\n']);
%     disp(error)
%     disp(indata.n_IR_l);
%     disp(max(error,[],'all'))
%     disp(k)
%     
% end


end