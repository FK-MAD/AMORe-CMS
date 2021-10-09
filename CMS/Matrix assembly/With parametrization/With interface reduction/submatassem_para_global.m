function submatdata=submatassem_para_global(submatdata,mats_S,mats,mats_interp_nom,dofdata,indata)
% submatdata -> variable in RAM OR matfile. works with both
% mats_S -> variable in RAM OR matfile. works with both
% mats -> variable in RAM OR matfile. works with both
% mats_interp_nom -> variable in RAM OR matfile. works with both

static=indata.static;

%% define additional needed quantities
T_tilde=dofdata.T_tilde;

n_I=dofdata.n_I;
n_i=dofdata.n_i;
n_id=indata.n_id;
n_IR=indata.n_IR;
N_S=indata.N_S;

S_0=indata.S_0;
S_j=indata.S_j;
func_g=indata.func_g;
func_h=indata.func_h;
n_theta=indata.n_theta;
theta_nom=indata.theta_nom;
invariant=indata.invariant;

%% without static correction


%% for substructures independent of theta

LAMBDA_id_S=mats_S.LAMBDA_id_S;
M_ib_S_hat=mats_S.M_ib_S_hat;

[M_ib_hat_0,LAMBDA_id_0]=diagassem(S_0,M_ib_S_hat,LAMBDA_id_S);

% (3.24) part of M_DI_hat independent of theta
submatdata.M_DI_hat_0=M_ib_hat_0*T_tilde;

% (3.25) part of K_DI_hat independent of theta
submatdata.K_DI_hat_0=LAMBDA_id_0;


%% for substructures that depend on theta

for k=1:n_theta    
    [M_ib_hat_j,LAMBDA_id_j]=diagassem(S_j{k},M_ib_S_hat,LAMBDA_id_S);
    
    % (3.24) part of M_DI_hat that depends on sqrt(g(theta))
    submatdata.M_DI_hat_j(1,k)={M_ib_hat_j*T_tilde};

    % (3.25) part of K_DI_hat that depends on h(theta)/g(theta)
    submatdata.K_DI_hat_j(1,k)={LAMBDA_id_j};
end


%% with static correction

if static==1
     
if invariant==1  % invariant assumption

F_S_bar=mats_S.F_S_bar;
F_bar=diagassem(S_0,F_S_bar);

LAMBDA_id=LAMBDA_id_0;
clear LAMBDA_id_0

% compute F_bar and LAMBDA_id at some nominal value of theta once
for k=1:n_theta
    [F_bar_j,LAMBDA_id_j]=diagassem(S_j{k},F_S_bar,LAMBDA_id_S);
    
    F_bar=F_bar+F_bar_j/func_h{k}(theta_nom(k));
    LAMBDA_id=LAMBDA_id+LAMBDA_id_j*func_h{k}(theta_nom(k))/func_g{k}(theta_nom(k));    
end
clear F_S_bar LAMBDA_id_S

YPSILON_I_theta_nom=mats_interp_nom.YPSILON_I_theta_k;

M_ib_hat=mats.M_ib_hat;

M_iIR=M_ib_hat*T_tilde*YPSILON_I_theta_nom;
clear M_ib_hat

M_ib_tilde=mats.M_ib_tilde;

L_bar=(F_bar*M_ib_tilde*T_tilde*YPSILON_I_theta_nom)/(speye(n_IR)-M_iIR'*M_iIR);
clear F_bar M_ib_tilde YPSILON_I_theta_nom

T_RI_1=-L_bar*M_iIR'*LAMBDA_id;
clear M_iIR LAMBDA_id

OMEGA_I_theta_nom=mats_interp_nom.OMEGA_I_theta_k;

T_RI_2=L_bar*OMEGA_I_theta_nom;
clear L_bar OMEGA_I_theta_nom

M_ii_S=mats_S.M_ii_S;
M_ib_S=mats_S.M_ib_S;
K_ii_S=mats_S.K_ii_S;
PSI_ib_S=mats_S.PSI_ib_S;
PHI_id_S=mats_S.PHI_id_S;

% structures holding needed expressions in (3.31) and (3.32)
temp1_S=cell(1,N_S);
temp2_S=cell(1,N_S);
temp3_S=cell(1,N_S);
for k=1:N_S
    temp1_S{k}=M_ii_S{k}*PHI_id_S{k};
    temp2_S{k}=M_ii_S{k}*PSI_ib_S{k}+M_ib_S{k};
    temp3_S{k}=K_ii_S{k}*PHI_id_S{k};
end
clear PHI_id_S PSI_ib_S M_ib_S M_ii_S

M_ii=mats.M_ii;
temp1=blkdiag(temp1_S{1:end-1},sparse(temp1_S{end}));
clear temp1_S
temp2=blkdiag(temp2_S{1:end-1},sparse(temp2_S{end}));
clear temp2_S

% submatrices of M_RI_hat

submatdata.M_RI_hat(1,1)={[T_RI_1'*M_ii*T_RI_1 T_RI_1'*M_ii*T_RI_2;...
    T_RI_2'*M_ii'*T_RI_1 T_RI_2'*M_ii*T_RI_2]};
clear M_ii

submatdata.M_RI_hat(1,2)={[T_RI_1'*temp1 sparse(n_id,n_IR);...
    T_RI_2'*temp1 sparse(n_IR,n_IR)]};
clear temp1

submatdata.M_RI_hat(1,3)={T_RI_1'*temp2*T_tilde};

submatdata.M_RI_hat(1,4)={T_RI_2'*temp2*T_tilde};
clear temp2

% submatrices of K_RI_hat

[K_ii_0,temp3_0]=diagassem(S_0,K_ii_S,temp3_S);

submatdata.K_RI_hat_0(1,1)={[T_RI_1'*K_ii_0*T_RI_1 T_RI_1'*K_ii_0*T_RI_2;...
    T_RI_2'*K_ii_0'*T_RI_1 T_RI_2'*K_ii_0*T_RI_2]};
clear K_ii_0

submatdata.K_RI_hat_0(1,2)={[T_RI_1'*temp3_0 sparse(n_id,n_IR);...
    T_RI_2'*temp3_0 sparse(n_IR,n_IR)]};    
clear temp3_0

for k=1:n_theta
    [K_ii_j,temp3_j]=diagassem(S_j{k},K_ii_S,temp3_S);
    
    submatdata.K_RI_hat_j(1,k)={[T_RI_1'*K_ii_j*T_RI_1 T_RI_1'*K_ii_j*T_RI_2;...
        T_RI_2'*K_ii_j'*T_RI_1 T_RI_2'*K_ii_j*T_RI_2]};
    clear K_ii_j
    
    submatdata.K_RI_hat_j(2,k)={[T_RI_1'*temp3_j sparse(n_id,n_IR);...
        T_RI_2'*temp3_j sparse(n_IR,n_IR)]};    
    clear temp3_j
end

    
else % full static correction
    

%% for substructures independent of theta
    
M_ii_S=mats_S.M_ii_S;
M_ib_S=mats_S.M_ib_S;
M_bb_S=mats_S.M_bb_S;
[M_ii_0,M_ib_0,M_bb_0]=diagassem(S_0,M_ii_S,M_ib_S,M_bb_S);
submatdata.M_hat_0=[M_ii_0 M_ib_0*T_tilde; T_tilde'*M_ib_0' T_tilde'*M_bb_0*T_tilde];
clear M_ii_0 M_ib_0 M_bb_0

K_ii_S=mats_S.K_ii_S;
K_ib_S=mats_S.K_ib_S;
K_bb_S=mats_S.K_bb_S;
[K_ii_0,K_ib_0,K_bb_0]=diagassem(S_0,K_ii_S,K_ib_S,K_bb_S);
submatdata.K_hat_0=[K_ii_0 K_ib_0*T_tilde; T_tilde'*K_ib_0' T_tilde'*K_bb_0*T_tilde];
clear K_ii_0 K_ib_0 K_bb_0

PHI_id_S=mats_S.PHI_id_S;
[PHI_id_0]=diagassem(S_0,PHI_id_S);

submatdata.T_DI_0(1,1)={PHI_id_0};
clear PHI_id_0

PSI_ib=mats.PSI_ib;

submatdata.T_DI_0(1,2)={PSI_ib*T_tilde};
clear PSI_ib

M_ib_S_tilde=mats_S.M_ib_S_tilde;
F_S_bar=mats_S.F_S_bar;
[submatdata.F_bar_0,submatdata.M_ib_tilde_0]=diagassem(S_0,F_S_bar,M_ib_S_tilde);

submatdata.M_iIR_0=M_ib_hat_0*T_tilde;
clear M_ib_hat_0
submatdata.LAMBDA_id_0=LAMBDA_id_0;
clear LAMBDA_id_0

%% for substructures that depend on theta

for k=1:n_theta    
    [M_ii_j,M_ib_j,M_bb_j]=diagassem(S_j{k},M_ii_S,M_ib_S,M_bb_S);
    submatdata.M_hat_j(1,k)={[M_ii_j M_ib_j*T_tilde; T_tilde'*M_ib_j' T_tilde'*M_bb_j*T_tilde]};
    clear M_ii_j M_ib_j M_bb_j
    
    [K_ii_j,K_ib_j,K_bb_j]=diagassem(S_j{k},K_ii_S,K_ib_S,K_bb_S);
    submatdata.K_hat_j(1,k)={[K_ii_j K_ib_j*T_tilde; T_tilde'*K_ib_j' T_tilde'*K_bb_j*T_tilde]};
    clear K_ii_j K_ib_j K_bb_j

    [PHI_id_j]=diagassem(S_j{k},PHI_id_S);
    submatdata.T_DI_j(1,k)={[PHI_id_j sparse(n_i,n_IR); sparse(n_I,n_id) sparse(n_I,n_IR)]};
    clear PHI_id_j
    
    [F_bar_j]=diagassem(S_j{k},F_S_bar);
    submatdata.F_bar_j(1,k)={F_bar_j};
    clear F_bar_j
    
    [M_ib_tilde_j]=diagassem(S_j{k},M_ib_S_tilde);
    submatdata.M_ib_tilde_j(1,k)={M_ib_tilde_j};
    clear M_ib_tilde_j

    [M_ib_hat_j,LAMBDA_id_j]=diagassem(S_j{k},M_ib_S_hat,LAMBDA_id_S);    
    submatdata.M_iIR_j(1,k)={M_ib_hat_j*T_tilde};
    submatdata.LAMBDA_id_j(1,k)={LAMBDA_id_j};
    clear M_ib_hat_j LAMBDA_id_j
end

end

end

end