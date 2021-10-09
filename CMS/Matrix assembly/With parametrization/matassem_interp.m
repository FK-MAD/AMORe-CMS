function mats_interp=matassem_interp(mats_interp,submatdata_interp,submatdata_M_I_K_I,dofdata,indata)
% mats_interp -> variable in RAM OR matfile. works with both
% submatdata_interp -> variable in RAM OR matfile. works with both
% submatdata_M_I_K_I -> variable in RAM OR matfile. works with both

theta_0=indata.theta_0;
theta_l=indata.theta_l; % a matrix of dimensions n_theta x L; each column contains the values of theta at a support point
theta_k=indata.theta_k;
L=indata.L; % # of support points

quad_interp=indata.quad_interp;

reduction_I=indata.reduction_I;

YPSILON_I_theta_0=submatdata_interp.YPSILON_I_theta_0;
YPSILON_I_theta_l=submatdata_interp.YPSILON_I_theta_l;
YPSILON_I_theta_ml=submatdata_interp.YPSILON_I_theta_ml;
theta_ml=submatdata_interp.theta_ml;

if reduction_I==5

    xi_k=find_weights(theta_0,theta_k,theta_l,theta_ml,L);

    YPSILON_Ill_l_hat_theta_k=assemble_modes_hat(xi_k,YPSILON_I_theta_0,YPSILON_I_theta_l,YPSILON_I_theta_ml,L,quad_interp,[]);
    
    [mats_interp.YPSILON_I_theta_k,mats_interp.OMEGA_I_theta_k]=assemble_modes(YPSILON_Ill_l_hat_theta_k,submatdata_M_I_K_I,dofdata,indata,[]);

elseif reduction_I==6
    
    I_dofs=dofdata.I_dofs;
    I_dofs_l=dofdata.I_dofs_l;
    N_I=dofdata.N_I;
    
    fprintf('\nInterpolating modes of all interfaces in parallel...\n');
    YPSILON_Ill_l_theta_k=cell(1,N_I);
    parfor k=1:N_I
        
        kept_param=find_kept_param(dofdata,indata,k);

        theta_0_temp=theta_0(kept_param);
        theta_l_temp=theta_l(kept_param,:);
        theta_k_temp=theta_k(kept_param);
        theta_ml_temp=[];
        if ~isempty(theta_ml)
            theta_ml_temp=theta_ml(kept_param,:);
        end

        if ~isempty(kept_param) % find kept parameters only if interface k depends at least on one parameter
            xi_k=find_weights(theta_0_temp,theta_k_temp,theta_l_temp,theta_ml_temp,L);
        else % if interface k does not depend on any parameters -> xi_k=[0,0,...,0]' to use only the nominal point which is exact
            xi_k=zeros(L,1);
        end
        
        YPSILON_Ill_l_hat_theta_k=assemble_modes_hat(xi_k,YPSILON_I_theta_0,YPSILON_I_theta_l,YPSILON_I_theta_ml,L,quad_interp,k);
             
        [YPSILON_Ill_l_theta_k{k},~]=assemble_modes(YPSILON_Ill_l_hat_theta_k,submatdata_M_I_K_I,dofdata,indata,k);

    end
    
    YPSILON_Ill_theta_k=blkdiag(YPSILON_Ill_l_theta_k{1:end-1},sparse(YPSILON_Ill_l_theta_k{end}));
    clear YPSILON_Ill_l_theta_k
    
    % re-sort rows of YPSILON_Ill to match numbering of interface dofs (in M_I and K_I) 
    [~,~,row_index]=intersect(I_dofs,vertcat(I_dofs_l{:}),'stable');
    mats_interp.YPSILON_I_theta_k=YPSILON_Ill_theta_k(row_index,:);
end

end

function xi_k=find_weights(theta_0,theta_k,theta_l,theta_ml,L)

    theta_l_ml=[theta_l theta_ml]; % matrix containing both theta(l) and theta(-l) (if exist)
    % size(theta_l_ml) ensures that if theta_ml is empty, repmat
    % creates a matrix of dimensions n_theta x L
    [~,nearest]=min(vecnorm(theta_l_ml-repmat(theta_k,1,size(theta_l_ml,2)),2,1));

    theta_q=theta_l_ml(:,nearest);

    xi_q_k=((theta_k-theta_0)'*(theta_q-theta_0))/norm(theta_q-theta_0)^2;

    % xi_q_k is NaN when theta_q=theta_0. If this happens, setting xi_q_k=1
    % results only in support point "nearest" to be acounted for.
    if isnan(xi_q_k)
        xi_q_k=1;
    end

    if nearest>L % if nearest point belongs to set theta(-l)
        xi_q_k=-xi_q_k; 
        nearest=nearest-L; % nearest must be between 1 and L when deleting the column of A that coresponds to theta_q
    end

    v_k=(theta_k-theta_0)-xi_q_k*(theta_q-theta_0);

    A=theta_l-repmat(theta_0,1,L);
    A(:,nearest)=[];

    tau_k=lsqminnorm(A,v_k);

    xi_k=[tau_k(1:nearest-1);xi_q_k;tau_k(nearest:end)];

% xi_k=lsqminnorm(A,v_k);
% xi_k(nearest)=xi_k(nearest)+xi_q_k;
end

function kept_param=find_kept_param(dofdata,indata,interface)
    adj_I=dofdata.adj_I;
    bound_S=dofdata.bound_S;
    group_l=indata.group_l;
    N_S=indata.N_S;
    S_0=indata.S_0;
    S_j=indata.S_j;
    boundaries=vertcat(adj_I{group_l{interface}});
    kept_param=[];
    for k=1:N_S 
        if ~isempty(intersect(boundaries,bound_S{k}))
            param=[];
            if ismember(k,S_0)
                %param=0;
                % if interface is adjacent to a group which does not depend to any model parameter,
                % do not add any parameter (zero has no meaning and breaks the code)
                param=[]; 
            else
                l=0;
                while isempty(param)
                    l=l+1;
                    if ismember(k,S_j{l})
                        param=l;
                    end
                end
            end
            kept_param=[kept_param param];
        end
    end
    kept_param=unique(kept_param);
end

function YPSILON_I_hat_theta_k=assemble_modes_hat(xi_k,YPSILON_I_theta_0,YPSILON_I_theta_l,YPSILON_I_theta_ml,L,quad_interp,interface)

if isempty(interface) % global with parametrization
    if quad_interp==0

        YPSILON_I_hat_theta_k=(1-sum(xi_k))*YPSILON_I_theta_0;

        for k=1:L
            YPSILON_I_hat_theta_k=YPSILON_I_hat_theta_k+xi_k(k)*YPSILON_I_theta_l{k};
        end

    else

        YPSILON_I_hat_theta_k=(1-sum(xi_k.^2))*YPSILON_I_theta_0;

        for k=1:L
            YPSILON_I_hat_theta_k=YPSILON_I_hat_theta_k+0.5*(xi_k(k)^2+xi_k(k))*YPSILON_I_theta_l{k}+0.5*(xi_k(k)^2-xi_k(k))*YPSILON_I_theta_ml{k};
        end

    end
else % local with parametrization
    if quad_interp==0
        
        YPSILON_I_hat_theta_k=(1-sum(xi_k))*YPSILON_I_theta_0{interface};
        
        for k=1:L
            YPSILON_I_hat_theta_k=YPSILON_I_hat_theta_k+xi_k(k)*YPSILON_I_theta_l{k}{interface};
        end
        
    else
        
         YPSILON_I_hat_theta_k=(1-sum(xi_k.^2))*YPSILON_I_theta_0{interface};

        for k=1:L
            YPSILON_I_hat_theta_k=YPSILON_I_hat_theta_k+0.5*(xi_k(k)^2+xi_k(k))*YPSILON_I_theta_l{k}{interface}+0.5*(xi_k(k)^2-xi_k(k))*YPSILON_I_theta_ml{k}{interface};
        end      
        
    end
       
end

end

function [YPSILON_I_theta_k,OMEGA_I_theta_k]=assemble_modes(YPSILON_I_hat_theta_k,submatdata_M_I_K_I,dofdata,indata,interface)

M_I_0=submatdata_M_I_K_I.M_I_0;
K_I_0=submatdata_M_I_K_I.K_I_0;
M_I_j=submatdata_M_I_K_I.M_I_j;
K_I_j=submatdata_M_I_K_I.K_I_j;

func_g=indata.func_g;
func_h=indata.func_h;
n_theta=indata.n_theta;
theta_k=indata.theta_k;

I_dofs_l=dofdata.I_dofs_l;
I_dofs=dofdata.I_dofs;

if ~isempty(interface) % local reduction
    [~,~,index]=intersect(I_dofs_l{interface},I_dofs,'stable');
    M_I_0=M_I_0(index,index);
    K_I_0=K_I_0(index,index);
    for k=1:n_theta
        M_I_j{k}=M_I_j{k}(index,index);
        K_I_j{k}=K_I_j{k}(index,index);  
    end
end

M_I_theta_k=M_I_0;
K_I_theta_k=K_I_0;

for k=1:n_theta
    M_I_theta_k=M_I_theta_k+M_I_j{k}*func_g{k}(theta_k(k));
    K_I_theta_k=K_I_theta_k+K_I_j{k}*func_h{k}(theta_k(k));    
end

M_I_reduced=YPSILON_I_hat_theta_k'*M_I_theta_k*YPSILON_I_hat_theta_k;
K_I_reduced=YPSILON_I_hat_theta_k'*K_I_theta_k*YPSILON_I_hat_theta_k;

% % for some reason, there is an error when M_I_reduced=1 (scalar). 
% % A fix is to remove sigma ('smallestabs')
% if M_I_reduced==1 
%     [Q_theta_k,OMEGA_I_theta_k]=eigs(K_I_reduced,M_I_reduced,length(K_I_reduced));
% else
%     [Q_theta_k,OMEGA_I_theta_k]=eigs(K_I_reduced,M_I_reduced,length(K_I_reduced),'smallestabs');
% end

% used eig instead of eigs because:
% 1) all modes need to be computed
% 2) matrices are usually small
% 3) there is no bug when M_I_reduced=1 (scalar) 
[Q_theta_k,OMEGA_I_theta_k]=eig(full(K_I_reduced),full(M_I_reduced));
OMEGA_I_theta_k=sparse(OMEGA_I_theta_k);

% Normalize Q_theta_k wrt the mass
for k=1:size(Q_theta_k,2)
    Q_theta_k(:,k)=Q_theta_k(:,k)/sqrt(Q_theta_k(:,k)'*M_I_reduced*Q_theta_k(:,k));
end

YPSILON_I_theta_k=YPSILON_I_hat_theta_k*Q_theta_k;


end