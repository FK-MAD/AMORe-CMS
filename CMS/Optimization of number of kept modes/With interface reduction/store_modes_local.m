function mats_I=store_modes_local(mats_I,mats,dofdata,indata)
% mats_I -> variable in RAM OR matfile. works with both
% mats -> variable in RAM OR matfile. works with both

T_tilde=dofdata.T_tilde;
I_dofs_l=dofdata.I_dofs_l;
N_I=dofdata.N_I;
I_dofs=dofdata.I_dofs;

M_bb_hat=mats.M_bb_hat;
M_I=T_tilde'*M_bb_hat*T_tilde;
mats_I.M_I=M_I;
clear M_bb_hat

K_bb_hat=mats.K_bb_hat;
K_I=T_tilde'*K_bb_hat*T_tilde;
mats_I.K_I=K_I;
clear K_bb_hat

n_IR_l=indata.n_IR_l_store;
method_eigf=indata.eigf.interface.method_store;
target_eigf=indata.eigf.interface.target_store;
max_eigf=indata.eigf.interface.max_store;
step_eigf=indata.eigf.interface.step_store;
init_eigf=indata.eigf.interface.init_store;

M_Ill_l_sliced=cell(1,N_I);
K_Ill_l_sliced=cell(1,N_I);
for k=1:N_I
    [~,~,index]=intersect(I_dofs_l{k},I_dofs,'stable');
    M_Ill_l_sliced{k}=M_I(index,index);
    K_Ill_l_sliced{k}=K_I(index,index);
end

    
fprintf('\nCalculating modes of all interfaces in parallel...\n');
parfor k=1:N_I
    [YPSILON_Ill_l,OMEGA_I,n_IR_l(k)]=keptmodes(K_Ill_l_sliced{k},M_Ill_l_sliced{k},method_eigf,n_IR_l(k),target_eigf(k),max_eigf(k),step_eigf(k),init_eigf(k));
    OMEGA_I_store{k}=sparse(OMEGA_I);
    OMEGA_I=[];
    
    % Normalize YPSILON_Ill_l wrt the mass
    for l=1:n_IR_l(k)
        YPSILON_Ill_l(:,l)=YPSILON_Ill_l(:,l)/sqrt(YPSILON_Ill_l(:,l)'*M_Ill_l_sliced{k}*YPSILON_Ill_l(:,l));
    end
    
    YPSILON_I_store{k}=YPSILON_Ill_l;
end
clear M_Ill_l_sliced K_Ill_l_sliced

mats_I.YPSILON_I_store=YPSILON_I_store;
mats_I.OMEGA_I_store=OMEGA_I_store;

end