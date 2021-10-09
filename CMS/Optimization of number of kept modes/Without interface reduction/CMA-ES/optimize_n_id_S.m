function optimize_n_id_S(dofdata,indata,LAMBDA)
% r=10.0186    7.7322    8.2230    9.6809    2.2649    3.6983
% 10.1449    6.0214 -> n_id_S=78    40    45     9     4     1    35     5

% N_S=indata.N_S;
% static=indata.static;
num_modes=indata.num_modes;
LAMBDA=LAMBDA(1:num_modes);

indata.eigf.group.method=1; % to keep modes based on a cutoff frequency

%mats=matfile('mats.mat','Writable',true);
mats=load('mats.mat');
% mats=[];
load('mats_S.mat','PSI_ib_S');
load('mats_S.mat','M_ii_S');
load('mats_S.mat','M_ib_S');
load('mats_S.mat','K_ii_S');
load('mats_S.mat','PHI_id_S_store');
load('mats_S.mat','LAMBDA_id_S_store');

tol=10^-3;

input_cma_opts
% Run CMAES to find the most probable values
[theta1,fval]=cmaes_v3_61('obj_n_id_S',theta0,[],opts,mats,PSI_ib_S,M_ii_S,M_ib_S,K_ii_S,PHI_id_S_store,LAMBDA_id_S_store,tol,dofdata,indata,LAMBDA);
  
% options = optimoptions('fmincon','algorithm','interior-point','OptimalityTolerance',1e-6,'StepTolerance',1e-10,...
%     'ConstraintTolerance',1e-6,'UseParallel',false,'Display','iter-detailed');

%  theta0=6*ones(indata.N_S,1);
%  lb=0.01*ones(indata.N_S,1);
%  ub=10*ones(indata.N_S,1);
% tol=10^-3;
% theta_check=[];
% error=[];

% [x,b]=fmincon(@obj_func,theta0,[],[],[],[],lb,ub,@cons,options);

% options = optimset('Display','final'); %iter
% 
% fminsearch(@obj_func_search,theta0,options);
% 
%     function f=obj_func_search(theta)
%         %theta=[5.5665    5.0787    5.3452    5.7122    6.0777    6.2530    6.0305    5.7130];
%         indata.eigf.group.target=theta'*2000;
%         
%         %tic
%         [mats,PHI_id_S_store,LAMBDA_id_S_store,indata]=...
%             update_mats_RAM(mats,PSI_ib_S,M_ii_S,M_ib_S,...
%             K_ii_S,PHI_id_S_store,LAMBDA_id_S_store,dofdata,indata);
%         %toc
%         
%         %tic
%         matdata=matassem_unreduced([],mats,dofdata,indata);
%         %toc
%         
%         %tic
%         LAMBDA_reduced=eigs(matdata.K_D_reduced,matdata.M_D_reduced,length(LAMBDA),'smallestabs');
%         %toc
%         
%         error=abs(sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)'    
%         
%         M=0;
%         if any(error>tol) || any(theta<lb) || any(theta>ub)
%             M=1000;
%         end
%         
%         disp(theta');
%         %disp(indata.n_id_S);
%         
%         f=sum(theta.^2)+M;
%     end

%     function f=obj_func(theta)
%         
%         f=sum(theta.^2);
%         disp(theta')
%         %disp(indata.n_id_S)
%         %f=sum(theta.^2)+sum(indata.n_id_S.^2);
%                         
%     end
% 
%     function [c,ceq]=cons(theta)
%         
%             indata.eigf.group.target=theta'*2000;
%             
%             %tic
%             [mats,PHI_id_S_store,LAMBDA_id_S_store,indata]=...
%                 update_mats_RAM(mats,PSI_ib_S,M_ii_S,M_ib_S,...
%                 K_ii_S,PHI_id_S_store,LAMBDA_id_S_store,dofdata,indata);
%             %toc
% 
%             %tic
%             matdata=matassem_unreduced([],mats,dofdata,indata);
%             %toc
% 
%             %tic
%             LAMBDA_reduced=eigs(matdata.K_D_reduced,matdata.M_D_reduced,length(LAMBDA),'smallestabs');
%             %toc
% 
%             error=max(abs(sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)'); 
%         
%         c=error-tol;
%         %disp(c);
%         ceq=[];
%     end
% l=0;
% for k=1:.5:10 % dummy loop to test speed
%     indata.eigf.group.target=k*ones(1,N_S)*4.6;
% 
% %     [indata]=update_mats_S_k(dofdata,indata);
% %     mats=update_mats(mats,indata);
%     
% %     [mats_S,indata]=update_mats_S_k_RAM(mats_S,dofdata,indata);
% %     mats=update_mats_RAM(mats,mats_S);
%     
%     tic
%     [mats,indata]=update_mats_HYBRID(mats,dofdata,indata);
%     toc
%     
%     tic
%     matdata=matassem_unreduced([],mats,dofdata,indata);
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
%     LAMBDA_reduced=eigs(K_reduced,M_reduced,num_modes,'smallestabs');
%         
%     error=abs((sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)');
%     
%     l=l+1;
%     fprintf(['\nIteration #',num2str(l),'\n']);
%     disp(error)
%     disp(indata.n_id_S);
%     disp(max(error,[],'all'))
% end


end