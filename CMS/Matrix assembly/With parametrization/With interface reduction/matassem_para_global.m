function matdata=matassem_para_global(matdata,submatdata,mats_interp_k,dofdata,indata)
% matdata -> variable in RAM OR matfile. works with both
% submatdata -> variable in RAM OR matfile. works with both
% mats_interp_k -> variable in RAM OR matfile. works with both

static=indata.static;
func_g=indata.func_g;
func_h=indata.func_h;
n_theta=indata.n_theta;
theta_k=indata.theta_k;
invariant=indata.invariant;

YPSILON_I_theta_k=mats_interp_k.YPSILON_I_theta_k;
OMEGA_I_theta_k=mats_interp_k.OMEGA_I_theta_k;

n_I=dofdata.n_I;
n_IR=indata.n_IR;
n_id=indata.n_id;

matdata.M_D_reduced=[speye(n_id) submatdata.M_DI_hat_0*YPSILON_I_theta_k;...
    YPSILON_I_theta_k'*submatdata.M_DI_hat_0' speye(n_IR)];

matdata.K_D_reduced=[submatdata.K_DI_hat_0 sparse(n_id,n_IR);...
    sparse(n_id,n_IR)' OMEGA_I_theta_k];

for k=1:n_theta
    matdata.M_D_reduced=matdata.M_D_reduced+...
        [sparse(n_id,n_id) cell2mat(submatdata.M_DI_hat_j(1,k))*YPSILON_I_theta_k;...
        YPSILON_I_theta_k'*cell2mat(submatdata.M_DI_hat_j(1,k))' sparse(n_IR,n_IR)]*sqrt(func_g{k}(theta_k(k)));

    matdata.K_D_reduced=matdata.K_D_reduced+...
        [cell2mat(submatdata.K_DI_hat_j(1,k)) sparse(n_id,n_IR);...
        sparse(n_id,n_IR)' sparse(n_IR,n_IR)]*func_h{k}(theta_k(k))/func_g{k}(theta_k(k));
end

% define in case there is no static correction
matdata.M_R_reduced=[];
matdata.K_R_reduced=[];

if static==1
    
if invariant==1 % invariant assumption
    
    % M_RI_hat
    
    matdata.M_R_reduced=cell2mat(submatdata.M_RI_hat(1,2))+...
        [sparse(n_id,n_id) cell2mat(submatdata.M_RI_hat(1,3))*YPSILON_I_theta_k;...
        sparse(n_IR,n_id) cell2mat(submatdata.M_RI_hat(1,4))*YPSILON_I_theta_k];
    
    matdata.M_R_reduced=matdata.M_R_reduced+matdata.M_R_reduced';
    
    matdata.M_R_reduced=matdata.M_R_reduced+cell2mat(submatdata.M_RI_hat(1,1));
    
    matdata.M_R_reduced=matdata.M_R_reduced+matdata.M_D_reduced;
    
    
    
    % K_RI_hat
    
    matdata.K_R_reduced=cell2mat(submatdata.K_RI_hat_0(1,2));
    
    for k=1:n_theta
        matdata.K_R_reduced=matdata.K_R_reduced+...
            cell2mat(submatdata.K_RI_hat_j(2,k))*func_h{k}(theta_k(k));              
    end
    
    matdata.K_R_reduced=matdata.K_R_reduced+matdata.K_R_reduced';
    
    matdata.K_R_reduced=matdata.K_R_reduced+cell2mat(submatdata.K_RI_hat_0(1,1));
    
    for k=1:n_theta
        matdata.K_R_reduced=matdata.K_R_reduced+...
            cell2mat(submatdata.K_RI_hat_j(1,k))*func_h{k}(theta_k(k));              
    end
    
    matdata.K_R_reduced=matdata.K_R_reduced+matdata.K_D_reduced;
        
else % full static correction
   
        T_tilde=dofdata.T_tilde;

        F_bar=submatdata.F_bar_0;
        M_ib_tilde=submatdata.M_ib_tilde_0;
        M_iIR=submatdata.M_iIR_0*YPSILON_I_theta_k;
        LAMBDA_id=submatdata.LAMBDA_id_0;
        
        for k=1:n_theta                    
            F_bar=F_bar+cell2mat(submatdata.F_bar_j(1,k))/func_h{k}(theta_k(k));
            M_ib_tilde=M_ib_tilde+cell2mat(submatdata.M_ib_tilde_j(1,k))*func_g{k}(theta_k(k));
            M_iIR=M_iIR+cell2mat(submatdata.M_iIR_j(1,k))*YPSILON_I_theta_k*sqrt(func_g{k}(theta_k(k)));
            LAMBDA_id=LAMBDA_id+cell2mat(submatdata.LAMBDA_id_j(1,k))*func_h{k}(theta_k(k))/func_g{k}(theta_k(k));
        end
        
        L_bar=(F_bar*M_ib_tilde*T_tilde*YPSILON_I_theta_k)/(speye(n_IR)-M_iIR'*M_iIR);
        clear F_bar M_ib_tilde
        
        T_RI=[-L_bar*M_iIR'*LAMBDA_id L_bar*OMEGA_I_theta_k; sparse(n_I,n_id) sparse(n_I,n_IR)];
        clear L_bar M_iIR LAMBDA_id OMEGA_I_theta_k
        
        M_hat=submatdata.M_hat_0;
        K_hat=submatdata.K_hat_0;   
        T_DI=[cell2mat(submatdata.T_DI_0(1,1)) cell2mat(submatdata.T_DI_0(1,2))*YPSILON_I_theta_k;...
            sparse(n_I,n_id) YPSILON_I_theta_k];
        clear YPSILON_I_theta_k
        
        for k=1:n_theta
            M_hat=M_hat+cell2mat(submatdata.M_hat_j(1,k))*func_g{k}(theta_k(k));
            K_hat=K_hat+cell2mat(submatdata.K_hat_j(1,k))*func_h{k}(theta_k(k));
            T_DI=T_DI+cell2mat(submatdata.T_DI_j(1,k))/sqrt(func_g{k}(theta_k(k)));
        end
        
        matdata.M_R_reduced=matdata.M_D_reduced+T_RI'*M_hat*T_DI+T_DI'*M_hat*T_RI+T_RI'*M_hat*T_RI;
        clear M_hat
        matdata.K_R_reduced=matdata.K_D_reduced+T_RI'*K_hat*T_DI+T_DI'*K_hat*T_RI+T_RI'*K_hat*T_RI;        
        clear K_hat T_RI T_DI
end


end


end