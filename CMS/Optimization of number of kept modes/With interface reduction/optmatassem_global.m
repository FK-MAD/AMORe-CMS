function [matdata,mats_I,indata]=optmatassem_global(matdata,mats,mats_I,dofdata,indata)
% matdata -> variable in RAM OR matfile. works with both
% mats -> variable in RAM OR matfile. works with both
% mats_I -> variable in RAM OR matfile. works with both

static=indata.static;

%% define additional needed quantities
n_I=dofdata.n_I;
T_tilde=dofdata.T_tilde;
n_id=indata.n_id;

%% without static correction

n_IR=indata.n_IR;
method_eigf=indata.eigf.interface.method;
target_eigf=indata.eigf.interface.target;
max_eigf=indata.eigf.interface.max;
step_eigf=indata.eigf.interface.step;
init_eigf=indata.eigf.interface.init;

YPSILON_I=mats_I.YPSILON_I_store;
OMEGA_I=mats_I.OMEGA_I_store;

% Get the needed modes from the stored ones. If they are not enough
% solve the eigenproblem
if method_eigf==0 && length(OMEGA_I)>=n_IR % contain at least as many modes as requested
    YPSILON_I=YPSILON_I(:,1:n_IR);
    OMEGA_I=OMEGA_I(1:n_IR,1:n_IR);
elseif method_eigf==1 && max(OMEGA_I,[],'all')>=(target_eigf(1)*2*pi)^2 % contain at least the modes corresponding to requested cutoff frequency
    kept=sum(diag(OMEGA_I)<(target_eigf(1)*2*pi)^2)+1;
    n_IR=kept;
    YPSILON_I=YPSILON_I(:,1:kept);
    OMEGA_I=OMEGA_I(1:kept,1:kept);
else % stored modes are not enough -> must calculate new eigenproblem
    warning('Not enough interface modes stored. Solving the eigenproblem...\n');

    M_I=mats_I.M_I;
    K_I=mats_I.K_I;
    
    % matrices YPSILON_I and OMEGA_I containing the kept interface modes
    fprintf('Calculating interface modes...\n');
    [YPSILON_I,OMEGA_I,n_IR]=keptmodes(K_I,M_I,method_eigf,n_IR,target_eigf(1),max_eigf(1),step_eigf(1),init_eigf(1));
    clear K_I
    OMEGA_I=sparse(OMEGA_I);

    % Normalize YPSILON_I wrt the mass
    for l=1:n_IR
        YPSILON_I(:,l)=YPSILON_I(:,l)/sqrt(YPSILON_I(:,l)'*M_I*YPSILON_I(:,l));
    end
    clear M_I
    
    % update storage since the calculated modes are more than those
    % stored
    mats_I.YPSILON_I_store=YPSILON_I;
    mats_I.OMEGA_I_store=OMEGA_I;
end

    
% update number of kept modes and dimensions of reduced matrices based on
% the kept modes -> update input
indata.n_IR=n_IR;
indata.n_DI=indata.n_id+indata.n_IR;
     
LAMBDA_id=mats.LAMBDA_id;

M_ib_hat=mats.M_ib_hat;

matdata.M_D_reduced=[speye(n_id) M_ib_hat*T_tilde*YPSILON_I;...
    YPSILON_I'*T_tilde'*M_ib_hat' speye(n_IR)];

matdata.K_D_reduced=[LAMBDA_id sparse(n_id,n_IR);...
    sparse(n_id,n_IR)' OMEGA_I];

%% with static correction

% define in case there is no static correction
matdata.M_R_reduced=[];
matdata.K_R_reduced=[];

if static==1
          
    M_ib_tilde=mats.M_ib_tilde;

    F_bar=mats.F_bar;
    
    M_iIR=M_ib_hat*T_tilde*YPSILON_I;
    clear M_ib_hat
    
    L_bar=(F_bar*M_ib_tilde*T_tilde*YPSILON_I)/(speye(n_IR)-M_iIR'*M_iIR);
    clear F_bar M_ib_tilde
     
    T_RI=[-L_bar*M_iIR'*LAMBDA_id L_bar*OMEGA_I; sparse(n_I,n_id) sparse(n_I,n_IR)];
    clear L_bar M_iIR LAMBDA_id OMEGA_I
    
    PHI_id=mats.PHI_id;

    PSI_ib=mats.PSI_ib;
    
    % transformation matrix that considers the effect of the dominant fixed-interface
    % normal modes and interface reduction
    T_DI=[PHI_id PSI_ib*T_tilde*YPSILON_I; sparse(n_I,n_id) YPSILON_I];
    clear YPSILON_I PHI_id PSI_ib
    
    M_ii=mats.M_ii;
    
    M_ib=mats.M_ib;
    
    M_bb=mats.M_bb;

    M_hat=[M_ii M_ib*T_tilde; T_tilde'*M_ib' T_tilde'*M_bb*T_tilde]; % sparse
    clear M_ii M_ib M_bb
    
    matdata.M_R_reduced=matdata.M_D_reduced+T_RI'*M_hat*T_DI+T_DI'*M_hat*T_RI+T_RI'*M_hat*T_RI;
    clear M_hat
    
    K_ii=mats.K_ii;
    
    K_ib=mats.K_ib;
    
    K_bb=mats.K_bb;
    
    K_hat=[K_ii K_ib*T_tilde; T_tilde'*K_ib' T_tilde'*K_bb*T_tilde]; % sparse
    clear K_ii K_ib K_bb  
    
    matdata.K_R_reduced=matdata.K_D_reduced+T_RI'*K_hat*T_DI+T_DI'*K_hat*T_RI+T_RI'*K_hat*T_RI;
    clear K_hat T_RI T_DI
end

end