function mats_I=store_modes_global(mats_I,mats,dofdata,indata)
% mats_I -> variable in RAM OR matfile. works with both
% mats -> variable in RAM OR matfile. works with both

T_tilde=dofdata.T_tilde;

M_bb_hat=mats.M_bb_hat;
M_I=T_tilde'*M_bb_hat*T_tilde;
mats_I.M_I=M_I;
clear M_bb_hat

K_bb_hat=mats.K_bb_hat;
K_I=T_tilde'*K_bb_hat*T_tilde;
mats_I.K_I=K_I;
clear K_bb_hat

n_IR=indata.n_IR_store;
method_eigf=indata.eigf.interface.method_store;
target_eigf=indata.eigf.interface.target_store;
max_eigf=indata.eigf.interface.max_store;
step_eigf=indata.eigf.interface.step_store;
init_eigf=indata.eigf.interface.init_store;

% matrices YPSILON_I and OMEGA_I containing the kept interface modes
fprintf('Calculating interface modes...\n');
[YPSILON_I,OMEGA_I,n_IR]=keptmodes(K_I,M_I,method_eigf,n_IR,target_eigf(1),max_eigf(1),step_eigf(1),init_eigf(1));
clear K_I
mats_I.OMEGA_I_store=sparse(OMEGA_I);

% Normalize YPSILON_I wrt the mass
for l=1:n_IR
    YPSILON_I(:,l)=YPSILON_I(:,l)/sqrt(YPSILON_I(:,l)'*M_I*YPSILON_I(:,l));
end
clear M_I

mats_I.YPSILON_I_store=YPSILON_I;
   
end