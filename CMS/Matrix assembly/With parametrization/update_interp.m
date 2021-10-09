function submatdata_interp_new=update_interp(submatdata_M_I_K_I,dofdata,indata)

func_g=indata.func_g;
func_h=indata.func_h;
n_theta=indata.n_theta;
theta_0=indata.theta_0;
theta_new=indata.theta_k; % the new support point is the sample point (which does not lie inside the convex hull of the existing support points)
quad_interp=indata.quad_interp;

% to ensure that the number of kept interface modes (global or local
% reduction) corresponds to the number of kept modes computed at all other
% (existing) support points (imported from input file submatdata_interp).
% Else, there is error in the assembly of reduced matrices from
% sub-matrices due to dimension mismatch.
indata.eigf.interface.method=0;

M_I_0=submatdata_M_I_K_I.M_I_0;
K_I_0=submatdata_M_I_K_I.K_I_0;
M_I_j=submatdata_M_I_K_I.M_I_j;
K_I_j=submatdata_M_I_K_I.K_I_j;

M_I_theta_l=M_I_0;
K_I_theta_l=K_I_0;

for l=1:n_theta
    M_I_theta_l=M_I_theta_l+M_I_j{l}*func_g{l}(theta_new(l));
    K_I_theta_l=K_I_theta_l+K_I_j{l}*func_h{l}(theta_new(l));
end

fprintf('\nCalculating interface modes at the new support point...\n');
[submatdata_interp_new.YPSILON_I_theta_new,~]=interface_modes(K_I_theta_l,M_I_theta_l,dofdata,indata);


% define in case there is no quadratic interpolation used
theta_new_ml=[];
YPSILON_I_theta_new_ml={};

if quad_interp==1
    theta_new_ml=theta_0-(theta_new-theta_0);
        
    M_I_theta_ml=M_I_0;
    K_I_theta_ml=K_I_0;
    
    for l=1:n_theta
        M_I_theta_ml=M_I_theta_ml+M_I_j{l}*func_g{l}(theta_new_ml(l));
        K_I_theta_ml=K_I_theta_ml+K_I_j{l}*func_h{l}(theta_new_ml(l));
    end
    
    fprintf('\nCalculating interface modes at the new minus support point (quadratic interpolation used)...\n');
    [YPSILON_I_theta_new_ml,~]=interface_modes(K_I_theta_ml,M_I_theta_ml,dofdata,indata);
end

submatdata_interp_new.theta_new_ml=theta_new_ml;
submatdata_interp_new.YPSILON_I_theta_new_ml=YPSILON_I_theta_new_ml;

end

function [YPSILON,indata]=interface_modes(K_I,M_I,dofdata,indata)

method_eigf=indata.eigf.interface.method;
target_eigf=indata.eigf.interface.target;
max_eigf=indata.eigf.interface.max;
step_eigf=indata.eigf.interface.step;
init_eigf=indata.eigf.interface.init;

reduction_I=indata.reduction_I;

if reduction_I==5 % global with parametrization
    n_IR=indata.n_IR;
    
    [YPSILON,~,n_IR]=keptmodes(K_I,M_I,method_eigf,n_IR,target_eigf(1),max_eigf(1),step_eigf(1),init_eigf(1));
    
    % Normalize YPSILON wrt the mass
    for k=1:size(YPSILON,2)
        YPSILON(:,k)=YPSILON(:,k)/sqrt(YPSILON(:,k)'*M_I*YPSILON(:,k));
    end

    indata.n_IR=n_IR;
    indata.n_DI=indata.n_id+indata.n_IR;
    
elseif reduction_I==6 % local with parametrization
    I_dofs_l=dofdata.I_dofs_l;
    I_dofs=dofdata.I_dofs;
    N_I=dofdata.N_I;
    n_IR_l=indata.n_IR_l;
    

    M_Ill_l_sliced=cell(1,N_I);
    K_Ill_l_sliced=cell(1,N_I);
    for k=1:N_I
        [~,~,index]=intersect(I_dofs_l{k},I_dofs,'stable');
        M_Ill_l_sliced{k}=M_I(index,index);
        K_Ill_l_sliced{k}=K_I(index,index);
    end
    
    fprintf('\nCalculating modes of all interfaces in parallel...\n');
    YPSILON=cell(1,N_I);
    parfor k=1:N_I
%         fprintf(['    Calculating modes of interface ',num2str(k),'...\n']);
%         [~,~,index]=intersect(I_dofs_l{k},I_dofs,'stable');
%         M_Ill_l=M_I(index,index);
%         K_Ill_l=K_I(index,index);
        [YPSILON{k},~,n_IR_l(k)]=keptmodes(K_Ill_l_sliced{k},M_Ill_l_sliced{k},method_eigf,n_IR_l(k),target_eigf(k),max_eigf(k),step_eigf(k),init_eigf(k));
        
        
        % Normalize YPSILON wrt the mass
        for l=1:size(YPSILON{k},2)
            YPSILON{k}(:,l)=YPSILON{k}(:,l)/sqrt(YPSILON{k}(:,l)'*M_Ill_l_sliced{k}*YPSILON{k}(:,l));
        end
        
    end
    clear M_Ill_l_sliced K_Ill_l_sliced
    
    % update number of kept modes and dimensions of reduced matrices based on
    % the kept modes -> update input
    indata.n_IR_l=n_IR_l;
    n_IRL=sum(n_IR_l);
    indata.n_IRL=n_IRL;
    indata.n_DIL=indata.n_id+indata.n_IRL;
    
end

end