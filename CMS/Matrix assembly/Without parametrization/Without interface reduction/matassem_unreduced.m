function matdata=matassem_unreduced(matdata,mats,dofdata,indata)
% matdata -> variable in RAM OR matfile. works with both
% mats -> variable in RAM OR matfile. works with both

static=indata.static;

%% define additional needed quantities
n_I=dofdata.n_I;
T_tilde=dofdata.T_tilde;
n_id=indata.n_id;

%% without static correction

M_ib_hat=mats.M_ib_hat;

M_bb_hat=mats.M_bb_hat;

K_bb_hat=mats.K_bb_hat;

LAMBDA_id=mats.LAMBDA_id;

matdata.M_D_reduced=[speye(n_id) M_ib_hat*T_tilde;...
    T_tilde'*M_ib_hat' T_tilde'*M_bb_hat*T_tilde]; % sparse

matdata.K_D_reduced=[LAMBDA_id sparse(n_id,n_I);...
    sparse(n_id,n_I)' T_tilde'*K_bb_hat*T_tilde]; % sparse

%% with static correction

% define in case there is no static correction
matdata.M_R_reduced=[];
matdata.K_R_reduced=[];

if static==1          
           
    M_ib_tilde=mats.M_ib_tilde;

    F_bar=mats.F_bar;
    
    M_I=T_tilde'*M_bb_hat*T_tilde;
    clear M_bb_hat
    M_iI=M_ib_hat*T_tilde;
    clear M_ib_hat
    
    L_bar=(F_bar*M_ib_tilde*T_tilde)/(M_I-M_iI'*M_iI);
    clear F_bar M_ib_tilde M_I
    
    K_I=T_tilde'*K_bb_hat*T_tilde;
    clear K_bb_hat
    
    T_R=[-L_bar*M_iI'*LAMBDA_id L_bar*K_I; sparse(n_I,n_id) sparse(n_I,n_I)];
    matdata.T_R_reduced=T_R;
    clear L_bar M_iI LAMBDA_id K_I
           
    PHI_id=mats.PHI_id;

    PSI_ib=mats.PSI_ib;

    T_D=[PHI_id PSI_ib*T_tilde; sparse(n_I,n_id) speye(n_I)]; % Craig-Bampton transformation matrix
    matdata.T_D_reduced=T_D;
    clear PHI_id PSI_ib

    M_ii=mats.M_ii;
    
    M_ib=mats.M_ib;
    
    M_bb=mats.M_bb;

    M_hat=[M_ii M_ib*T_tilde; T_tilde'*M_ib' T_tilde'*M_bb*T_tilde]; % sparse
    clear M_ii M_ib M_bb
       
    matdata.M_R_reduced=matdata.M_D_reduced+T_R'*M_hat*T_D+T_D'*M_hat*T_R+T_R'*M_hat*T_R;
    clear M_hat
    
    K_ii=mats.K_ii;
    
    K_ib=mats.K_ib;
    
    K_bb=mats.K_bb;
    
    K_hat=[K_ii K_ib*T_tilde; T_tilde'*K_ib' T_tilde'*K_bb*T_tilde]; % sparse
    clear K_ii K_ib K_bb
    
    matdata.K_R_reduced=matdata.K_D_reduced+T_R'*K_hat*T_D+T_D'*K_hat*T_R+T_R'*K_hat*T_R;    
    clear K_hat T_R T_D
end

end