function matdata=matassem_para_local(matdata,submatdata,mats_interp_k,dofdata,indata)
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

n_IRL=indata.n_IRL;
n_id=indata.n_id;

matdata.M_D_reduced=[speye(n_id) cell2mat(submatdata.M_DIL_hat_0(1,1))*YPSILON_I_theta_k;...
    YPSILON_I_theta_k'*cell2mat(submatdata.M_DIL_hat_0(1,1))' YPSILON_I_theta_k'*cell2mat(submatdata.M_DIL_hat_0(1,2))*YPSILON_I_theta_k];

matdata.K_D_reduced=[cell2mat(submatdata.K_DIL_hat_0(1,1)) sparse(n_id,n_IRL);...
    sparse(n_id,n_IRL)' YPSILON_I_theta_k'*cell2mat(submatdata.K_DIL_hat_0(1,2))*YPSILON_I_theta_k];

for k=1:n_theta
    matdata.M_D_reduced=matdata.M_D_reduced+...       
        [sparse(n_id,n_id) cell2mat(submatdata.M_DIL_hat_j(1,k))*YPSILON_I_theta_k;...
        YPSILON_I_theta_k'*cell2mat(submatdata.M_DIL_hat_j(1,k))' sparse(n_IRL,n_IRL)]*sqrt(func_g{k}(theta_k(k)))+...
        [sparse(n_id,n_id) sparse(n_id,n_IRL);...
        sparse(n_id,n_IRL)' YPSILON_I_theta_k'*cell2mat(submatdata.M_DIL_hat_j(2,k))*YPSILON_I_theta_k]*func_g{k}(theta_k(k));

    matdata.K_D_reduced=matdata.K_D_reduced+...
        [cell2mat(submatdata.K_DIL_hat_j(1,k)) sparse(n_id,n_IRL);...
        sparse(n_id,n_IRL)' sparse(n_IRL,n_IRL)]*func_h{k}(theta_k(k))/func_g{k}(theta_k(k))+...
        [sparse(n_id,n_id) sparse(n_id,n_IRL);...
        sparse(n_id,n_IRL)' YPSILON_I_theta_k'*cell2mat(submatdata.K_DIL_hat_j(2,k))*YPSILON_I_theta_k]*func_h{k}(theta_k(k));
end

% define in case there is no static correction
matdata.M_R_reduced=[];
matdata.K_R_reduced=[];

if static==1
    
    if invariant==1 % invariant assumption
        
        % M_RIL_hat
        
        matdata.M_R_reduced=cell2mat(submatdata.M_RIL_hat(1,2))+...
            [sparse(n_id,n_id) cell2mat(submatdata.M_RIL_hat(1,3))*YPSILON_I_theta_k;...
            sparse(n_IRL,n_id) cell2mat(submatdata.M_RIL_hat(1,4))*YPSILON_I_theta_k];
        
        matdata.M_R_reduced=matdata.M_R_reduced+matdata.M_R_reduced';
        
        matdata.M_R_reduced=matdata.M_R_reduced+cell2mat(submatdata.M_RIL_hat(1,1));
        
        matdata.M_R_reduced=matdata.M_R_reduced+matdata.M_D_reduced;
        
        
        
        % K_RIL_hat
        
        matdata.K_R_reduced=cell2mat(submatdata.K_RIL_hat_0(1,2));
        
        for k=1:n_theta
            matdata.K_R_reduced=matdata.K_R_reduced+...
                cell2mat(submatdata.K_RIL_hat_j(2,k))*func_h{k}(theta_k(k));
        end
        
        matdata.K_R_reduced=matdata.K_R_reduced+matdata.K_R_reduced';
        
        matdata.K_R_reduced=matdata.K_R_reduced+cell2mat(submatdata.K_RIL_hat_0(1,1));
        
        for k=1:n_theta
            matdata.K_R_reduced=matdata.K_R_reduced+...
                cell2mat(submatdata.K_RIL_hat_j(1,k))*func_h{k}(theta_k(k));
        end
        
        matdata.K_R_reduced=matdata.K_R_reduced+matdata.K_D_reduced;
        
    else % full static correction
        
        n_I=dofdata.n_I;
        T_tilde=dofdata.T_tilde;
        
        F_bar=submatdata.F_bar_0;
        M_ib_tilde=submatdata.M_ib_tilde_0;
        M_II=YPSILON_I_theta_k'*submatdata.M_II_0*YPSILON_I_theta_k;
        M_iIRL=submatdata.M_iIRL_0*YPSILON_I_theta_k;
        
        for k=1:n_theta
            F_bar=F_bar+cell2mat(submatdata.F_bar_j(1,k))/func_h{k}(theta_k(k));
            M_ib_tilde=M_ib_tilde+cell2mat(submatdata.M_ib_tilde_j(1,k))*func_g{k}(theta_k(k));
            M_II=M_II+YPSILON_I_theta_k'*cell2mat(submatdata.M_II_j(1,k))*YPSILON_I_theta_k*func_g{k}(theta_k(k));
            M_iIRL=M_iIRL+cell2mat(submatdata.M_iIRL_j(1,k))*YPSILON_I_theta_k*sqrt(func_g{k}(theta_k(k)));
        end
        
        L_bar=(F_bar*M_ib_tilde*T_tilde*YPSILON_I_theta_k)/(M_II-M_iIRL'*M_iIRL);
        clear F_bar M_ib_tilde M_II
        
        
        LAMBDA_id=submatdata.LAMBDA_id_0;
        K_II=YPSILON_I_theta_k'*submatdata.K_II_0*YPSILON_I_theta_k;
        for k=1:n_theta
            LAMBDA_id=LAMBDA_id+cell2mat(submatdata.LAMBDA_id_j(1,k))*func_h{k}(theta_k(k))/func_g{k}(theta_k(k));
            K_II=K_II+YPSILON_I_theta_k'*cell2mat(submatdata.K_II_j(1,k))*YPSILON_I_theta_k*func_h{k}(theta_k(k));
        end
        
        T_RIL=[-L_bar*M_iIRL'*LAMBDA_id L_bar*K_II; sparse(n_I,n_id) sparse(n_I,n_IRL)];
        clear L_bar M_iIRL LAMBDA_id K_II
        
        M_hat=submatdata.M_hat_0;
        K_hat=submatdata.K_hat_0;
        T_DIL=[cell2mat(submatdata.T_DIL_0(1,1)) cell2mat(submatdata.T_DIL_0(1,2))*YPSILON_I_theta_k;...
            sparse(n_I,n_id) YPSILON_I_theta_k];
        
        for k=1:n_theta
            M_hat=M_hat+cell2mat(submatdata.M_hat_j(1,k))*func_g{k}(theta_k(k));
            K_hat=K_hat+cell2mat(submatdata.K_hat_j(1,k))*func_h{k}(theta_k(k));
            T_DIL=T_DIL+cell2mat(submatdata.T_DIL_j(1,k))/sqrt(func_g{k}(theta_k(k)));
        end
        
        matdata.M_R_reduced=matdata.M_D_reduced+T_RIL'*M_hat*T_DIL+T_DIL'*M_hat*T_RIL+T_RIL'*M_hat*T_RIL;
        clear M_hat
        matdata.K_R_reduced=matdata.K_D_reduced+T_RIL'*K_hat*T_DIL+T_DIL'*K_hat*T_RIL+T_RIL'*K_hat*T_RIL;
        clear K_hat T_RIL T_DIL
        
    end

end

end