function matdata=matassem_para_unreduced(matdata,submatdata,dofdata,indata)
% matdata -> variable in RAM OR matfile. works with both
% submatdata -> variable in RAM OR matfile. works with both

static=indata.static;
func_g=indata.func_g;
func_h=indata.func_h;
n_theta=indata.n_theta;
theta_k=indata.theta_k;
invariant=indata.invariant;

matdata.M_D_reduced=submatdata.M_D_hat_0;
matdata.K_D_reduced=submatdata.K_D_hat_0;

for k=1:n_theta
    matdata.M_D_reduced=matdata.M_D_reduced+cell2mat(submatdata.M_D_hat_j(1,k))*sqrt(func_g{k}(theta_k(k)))+...
        cell2mat(submatdata.M_D_hat_j(2,k))*func_g{k}(theta_k(k));
    
    matdata.K_D_reduced=matdata.K_D_reduced+cell2mat(submatdata.K_D_hat_j(1,k))*func_h{k}(theta_k(k))/func_g{k}(theta_k(k))+...
        cell2mat(submatdata.K_D_hat_j(2,k))*func_h{k}(theta_k(k));
end

% define in case there is no static correction
matdata.M_R_reduced=[];
matdata.K_R_reduced=[];

if static==1
    
    n_I=dofdata.n_I;
    n_id=indata.n_id;
    
    if invariant==1 % invariant assumption
        
       T_R1=submatdata.T_R1;
       T_R2=submatdata.T_R2_0;
       
       for k=1:n_theta
            T_R2=T_R2+cell2mat(submatdata.T_R2_j(1,k))*func_h{k}(theta_k(k));
       end
       
       % M_R_hat
       
       matdata.M_R_reduced=[T_R1'*cell2mat(submatdata.M_R_hat(1,2)) sparse(n_id,n_I);...
           T_R2'*cell2mat(submatdata.M_R_hat(1,2)) sparse(n_I,n_I)]+...
           [sparse(n_id,n_id) T_R1'*cell2mat(submatdata.M_R_hat(1,3));...
           sparse(n_I,n_id) T_R2'*cell2mat(submatdata.M_R_hat(1,3))];
       
       matdata.M_R_reduced=matdata.M_R_reduced+matdata.M_R_reduced';
       
       matdata.M_R_reduced=matdata.M_R_reduced+...
           [T_R1'*cell2mat(submatdata.M_R_hat(1,1))*T_R1 T_R1'*cell2mat(submatdata.M_R_hat(1,1))*T_R2;...
           T_R2'*cell2mat(submatdata.M_R_hat(1,1))'*T_R1 T_R2'*cell2mat(submatdata.M_R_hat(1,1))*T_R2];
       
       matdata.M_R_reduced=matdata.M_R_reduced+matdata.M_D_reduced;
       
       
       
       % K_R_hat
       
       matdata.K_R_reduced=[T_R1'*cell2mat(submatdata.K_R_hat_0(1,2)) sparse(n_id,n_I);...
           T_R2'*cell2mat(submatdata.K_R_hat_0(1,2)) sparse(n_I,n_I)];
       
       for k=1:n_theta
           matdata.K_R_reduced=matdata.K_R_reduced+...
               [T_R1'*cell2mat(submatdata.K_R_hat_j(2,k)) sparse(n_id,n_I);...
               T_R2'*cell2mat(submatdata.K_R_hat_j(2,k)) sparse(n_I,n_I)]*func_h{k}(theta_k(k));       
       end
       
       matdata.K_R_reduced=matdata.K_R_reduced+matdata.K_R_reduced';
       
       matdata.K_R_reduced=matdata.K_R_reduced+...
           [T_R1'*cell2mat(submatdata.K_R_hat_0(1,1))*T_R1 T_R1'*cell2mat(submatdata.K_R_hat_0(1,1))*T_R2;...
           T_R2'*cell2mat(submatdata.K_R_hat_0(1,1))'*T_R1 T_R2'*cell2mat(submatdata.K_R_hat_0(1,1))*T_R2];
       
       for k=1:n_theta
           matdata.K_R_reduced=matdata.K_R_reduced+...
               [T_R1'*cell2mat(submatdata.K_R_hat_j(1,k))*T_R1 T_R1'*cell2mat(submatdata.K_R_hat_j(1,k))*T_R2;...
               T_R2'*cell2mat(submatdata.K_R_hat_j(1,k))'*T_R1 T_R2'*cell2mat(submatdata.K_R_hat_j(1,k))*T_R2]*func_h{k}(theta_k(k));
           
       end
       
       matdata.K_R_reduced=matdata.K_R_reduced+matdata.K_D_reduced;
       
    else % full static correction
        
        T_tilde=dofdata.T_tilde;
           
        F_bar=submatdata.F_bar_0;
        M_ib_tilde=submatdata.M_ib_tilde_0;
        M_I=submatdata.M_I_0;
        M_iI=submatdata.M_iI_0;

        for k=1:n_theta           
            F_bar=F_bar+cell2mat(submatdata.F_bar_j(1,k))/func_h{k}(theta_k(k));
            M_ib_tilde=M_ib_tilde+cell2mat(submatdata.M_ib_tilde_j(1,k))*func_g{k}(theta_k(k));
            M_I=M_I+cell2mat(submatdata.M_I_j(1,k))*func_g{k}(theta_k(k));
            M_iI=M_iI+cell2mat(submatdata.M_iI_j(1,k))*sqrt(func_g{k}(theta_k(k)));
        end
        
        L_bar=(F_bar*M_ib_tilde*T_tilde)/(M_I-M_iI'*M_iI);
        clear F_bar M_ib_tilde M_I
        
        LAMBDA_id=submatdata.LAMBDA_id_0;
        K_I=submatdata.K_I_0;
        for k=1:n_theta
            LAMBDA_id=LAMBDA_id+cell2mat(submatdata.LAMBDA_id_j(1,k))*func_h{k}(theta_k(k))/func_g{k}(theta_k(k));
            K_I=K_I+cell2mat(submatdata.K_I_j(1,k))*func_h{k}(theta_k(k)); 
        end
        
        T_R=[-L_bar*M_iI'*LAMBDA_id L_bar*K_I; sparse(n_I,n_id) sparse(n_I,n_I)];
        clear L_bar M_iI LAMBDA_id K_I
        
        M_hat=submatdata.M_hat_0;
        K_hat=submatdata.K_hat_0;
        T_D=submatdata.T_D_0;
        
        for k=1:n_theta
            M_hat=M_hat+cell2mat(submatdata.M_hat_j(1,k))*func_g{k}(theta_k(k));
            K_hat=K_hat+cell2mat(submatdata.K_hat_j(1,k))*func_h{k}(theta_k(k));
            T_D=T_D+cell2mat(submatdata.T_D_j(1,k))/sqrt(func_g{k}(theta_k(k)));   
        end
        
        matdata.M_R_reduced=matdata.M_D_reduced+T_R'*M_hat*T_D+T_D'*M_hat*T_R+T_R'*M_hat*T_R;
        clear M_hat
        matdata.K_R_reduced=matdata.K_D_reduced+T_R'*K_hat*T_D+T_D'*K_hat*T_R+T_R'*K_hat*T_R;
        clear K_hat T_R T_D
    end
       
end

end