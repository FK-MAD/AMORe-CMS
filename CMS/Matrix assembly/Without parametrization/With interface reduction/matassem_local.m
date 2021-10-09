function [matdata,indata]=matassem_local(matdata,mats,dofdata,indata)
% matdata -> variable in RAM OR matfile. works with both
% mats -> variable in RAM OR matfile. works with both

static=indata.static;

%% define additional needed quantities

n_I=dofdata.n_I;
T_tilde=dofdata.T_tilde;
I_dofs_l=dofdata.I_dofs_l;
N_I=dofdata.N_I;
I_dofs=dofdata.I_dofs;
n_id=indata.n_id;

%% without static correction

M_bb_hat=mats.M_bb_hat;

M_I=T_tilde'*M_bb_hat*T_tilde;
clear M_bb_hat

K_bb_hat=mats.K_bb_hat;

K_I=T_tilde'*K_bb_hat*T_tilde;
clear K_bb_hat

n_IR_l=indata.n_IR_l;
method_eigf=indata.eigf.interface.method;
target_eigf=indata.eigf.interface.target;
max_eigf=indata.eigf.interface.max;
step_eigf=indata.eigf.interface.step;
init_eigf=indata.eigf.interface.init;


M_Ill_l_sliced=cell(1,N_I);
K_Ill_l_sliced=cell(1,N_I);
for k=1:N_I
    [~,~,index]=intersect(I_dofs_l{k},I_dofs,'stable');
    M_Ill_l_sliced{k}=M_I(index,index);
    K_Ill_l_sliced{k}=K_I(index,index);
end

fprintf('\nCalculating modes of all interfaces in parallel...\n');
YPSILON_Ill_l=cell(1,N_I);
parfor k=1:N_I
    %fprintf(['Calculating modes of interface ',num2str(k),'...\n']);
    [YPSILON_Ill_l{k},~,n_IR_l(k)]=keptmodes(K_Ill_l_sliced{k},M_Ill_l_sliced{k},method_eigf,n_IR_l(k),target_eigf(k),max_eigf(k),step_eigf(k),init_eigf(k));

    % Normalize YPSILON_Ill_l wrt the mass
    for l=1:n_IR_l(k)
        YPSILON_Ill_l{k}(:,l)=YPSILON_Ill_l{k}(:,l)/sqrt(YPSILON_Ill_l{k}(:,l)'*M_Ill_l_sliced{k}*YPSILON_Ill_l{k}(:,l));
    end
end
clear M_Ill_l_sliced K_Ill_l_sliced

% update number of kept modes and dimensions of reduced matrices based on
% the kept modes -> update input
indata.n_IR_l=n_IR_l;
n_IRL=sum(n_IR_l);
indata.n_IRL=n_IRL;
indata.n_DIL=indata.n_id+indata.n_IRL;

YPSILON_Ill=blkdiag(YPSILON_Ill_l{1:end-1},sparse(YPSILON_Ill_l{end})); % this way blkdiag returns sparse matrix 
clear YPSILON_Ill_l

% %% IS IT NEEDED? -> NO!
% % re-sort rows of YPSILON_Ill to match numbering of interface dofs (in M_I and K_I) 
% [~,~,row_index]=intersect(I_dofs,vertcat(I_dofs_l{:}),'stable');
% YPSILON_Ill=YPSILON_Ill(row_index,:);
% %% 

% e.g.
% YPSILON_Ill=
%      1     1     0     0
%      2     2     0     0
%      3     3     0     0
%      6     6     0     0
%      7     7     0     0
%      0     0     4     4
%      0     0     5     5
%      0     0     8     8
%      0     0     9     9
%      0     0    10    10
% becomes
% YPSILON_Ill=
%      1     1     0     0
%      2     2     0     0
%      3     3     0     0
%      0     0     4     4
%      0     0     5     5
%      6     6     0     0
%      7     7     0     0
%      0     0     8     8
%      0     0     9     9
%      0     0    10    10

K_II=YPSILON_Ill'*K_I*YPSILON_Ill;
clear K_I
M_II=YPSILON_Ill'*M_I*YPSILON_Ill;
clear M_I

M_ib_hat=mats.M_ib_hat;

M_iIRL=M_ib_hat*T_tilde*YPSILON_Ill;
clear M_ib_hat
    
LAMBDA_id=mats.LAMBDA_id;

matdata.M_D_reduced=[speye(n_id) M_iIRL; M_iIRL' M_II];
matdata.K_D_reduced=[LAMBDA_id sparse(n_id,n_IRL); sparse(n_id,n_IRL)' K_II];

%% with static correction

% define in case there is no static correction
matdata.M_R_reduced=[];
matdata.K_R_reduced=[];

if static==1
       
    M_ib_tilde=mats.M_ib_tilde;

    F_bar=mats.F_bar;
    
    L_bar=(F_bar*M_ib_tilde*T_tilde*YPSILON_Ill)/(M_II-M_iIRL'*M_iIRL);
    clear F_bar M_ib_tilde M_II
    
    T_RIL=[-L_bar*M_iIRL'*LAMBDA_id L_bar*K_II; sparse(n_I,n_id) sparse(n_I,n_IRL)];
    matdata.T_R_reduced=T_RIL;
    clear L_bar M_iIRL LAMBDA_id K_II
    
    PHI_id=mats.PHI_id;

    PSI_ib=mats.PSI_ib;
        
    % transformation matrix that considers the effect of the dominant fixed-interface
    % normal modes and local interface reduction
    T_DIL=[PHI_id PSI_ib*T_tilde*YPSILON_Ill; sparse(n_I,n_id) YPSILON_Ill];
    matdata.T_D_reduced=T_DIL;
    clear YPSILON_Ill PHI_id PSI_ib
    
    M_ii=mats.M_ii;
    
    M_ib=mats.M_ib;
    
    M_bb=mats.M_bb;

    M_hat=[M_ii M_ib*T_tilde; T_tilde'*M_ib' T_tilde'*M_bb*T_tilde]; % sparse
    clear M_ii M_ib M_bb
    
    matdata.M_R_reduced=matdata.M_D_reduced+T_RIL'*M_hat*T_DIL+T_DIL'*M_hat*T_RIL+T_RIL'*M_hat*T_RIL;
    clear M_hat
    
    K_ii=mats.K_ii;
    
    K_ib=mats.K_ib;
    
    K_bb=mats.K_bb;
    
    K_hat=[K_ii K_ib*T_tilde; T_tilde'*K_ib' T_tilde'*K_bb*T_tilde]; % sparse
    clear K_ii K_ib K_bb   
    
    matdata.K_R_reduced=matdata.K_D_reduced+T_RIL'*K_hat*T_DIL+T_DIL'*K_hat*T_RIL+T_RIL'*K_hat*T_RIL;    
    clear K_hat T_RIL T_DIL
end

end