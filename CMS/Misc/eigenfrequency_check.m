Kc=full(moddata.mats.Kc);
Mc=full(moddata.mats.Ec);
lambda_c=eigs(Kc,Mc,6,'smallestabs');
sqrt(lambda_c)'

K=full(moddata.mats.K);
M=full(moddata.mats.E);
lambda=eigs(K,M,10,'smallestabs');
sqrt(lambda)'

K_hat=full(matdata.K_hat);
M_hat=full(matdata.M_hat);
lambda_hat=eigs(K_hat,M_hat,6,'smallestabs');
sqrt(lambda_hat)'

T_G=full(dofdata.T_G);
M_hat_perm=T_G*M_hat*T_G';
K_hat_perm=T_G*K_hat*T_G';
lambda_hat_perm=eigs(K_hat_perm,M_hat_perm,10,'smallestabs');
sqrt(lambda_hat_perm)'

K_D=full(matdata.K_D_reduced);
M_D=full(matdata.M_D_reduced);
K_D=(K_D'+K_D)/2;
M_D=(M_D'+M_D)/2;
lambda_D=eigs(K_D,M_D,6,'smallestabs');
sqrt(lambda_D)'

K_R=full(matdata.K_R_reduced);
M_R=full(matdata.M_R_reduced);
lambda_R=eigs(K_R,M_R,6,'smallestabs');
sqrt(lambda_R)'

M_D_perm=T_G*M_D*T_G';
K_D_perm=T_G*K_D*T_G';
lambda_D_perm=eigs(K_D_perm,M_D_perm,6,'smallestabs');
sqrt(lambda_D_perm)'
