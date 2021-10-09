function submatdata_M_I_K_I=submatassem_M_I_K_I(submatdata_M_I_K_I,mats_S,dofdata,indata)
% submatdata_M_I_K_I -> variable in RAM OR matfile. works with both

K_bb_S_hat=mats_S.K_bb_S_hat;
M_bb_S_hat=mats_S.M_bb_S_hat;

T_tilde=dofdata.T_tilde;

S_0=indata.S_0;
S_j=indata.S_j;

n_theta=indata.n_theta;

[M_bb_hat_0,K_bb_hat_0]=diagassem(S_0,M_bb_S_hat,K_bb_S_hat);

submatdata_M_I_K_I.M_I_0=T_tilde'*M_bb_hat_0*T_tilde;
submatdata_M_I_K_I.K_I_0=T_tilde'*K_bb_hat_0*T_tilde;

for k=1:n_theta    
    [M_bb_hat_j,K_bb_hat_j]=diagassem(S_j{k},M_bb_S_hat,K_bb_S_hat);
    submatdata_M_I_K_I.M_I_j(1,k)={T_tilde'*M_bb_hat_j*T_tilde};
    submatdata_M_I_K_I.K_I_j(1,k)={T_tilde'*K_bb_hat_j*T_tilde};    
end


end