function [submatdata_interp,indata]=submatassem_interp(submatdata_interp,submatdata_M_I_K_I,dofdata,indata)
% submatdata_interp -> variable in RAM OR matfile. works with both
% submatdata_M_I_K_I -> variable in RAM OR matfile. works with both

func_g=indata.func_g;
func_h=indata.func_h;
n_theta=indata.n_theta;
theta_0=indata.theta_0;
theta_l=indata.theta_l; % a matrix of dimensions n_theta x L; each column contains the values of theta at a support point
quad_interp=indata.quad_interp;
L=indata.L; % # of support points

M_I_0=submatdata_M_I_K_I.M_I_0;
K_I_0=submatdata_M_I_K_I.K_I_0;
M_I_j=submatdata_M_I_K_I.M_I_j;
K_I_j=submatdata_M_I_K_I.K_I_j;

M_I_theta_0=M_I_0;
K_I_theta_0=K_I_0;

for k=1:n_theta
    M_I_theta_0=M_I_theta_0+M_I_j{k}*func_g{k}(theta_0(k));
    K_I_theta_0=K_I_theta_0+K_I_j{k}*func_h{k}(theta_0(k));    
end

fprintf('Calculating interface modes at nominal point...\n');
[submatdata_interp.YPSILON_I_theta_0,indata]=interface_modes(K_I_theta_0,M_I_theta_0,dofdata,indata);
clear M_I_theta_0 K_I_theta_0

% ensures that for each support point the same number of modes is kept as
% with the nominal point. Otherwise there is a mismatch of dimensions in
% the assembly.
indata.eigf.interface.method=0;


% create a waitbar to observe progress
D = parallel.pool.DataQueue;
h = waitbar(0,'Progress of parallel computations','Name','Please wait ...');
p=1;
afterEach(D, @nUpdateWaitbar);


fprintf('\nCalculating interface modes at all support points in parallel...\n');
YPSILON_I_theta_l=cell(1,L);
parfor k=1:L  
    M_I_theta_l=M_I_0;
    K_I_theta_l=K_I_0;
    
    for l=1:n_theta
        M_I_theta_l=M_I_theta_l+M_I_j{l}*func_g{l}(theta_l(l,k));
        K_I_theta_l=K_I_theta_l+K_I_j{l}*func_h{l}(theta_l(l,k));
    end
    
    fprintf(['\nCalculating interface modes at support point ',num2str(k),'...\n']);
    [YPSILON_I_theta_l{k},~]=interface_modes(K_I_theta_l,M_I_theta_l,dofdata,indata);
    
    send(D, k); % update waitbar
end
delete(h);
clear K_I_theta_l M_I_theta_l


submatdata_interp.YPSILON_I_theta_l=YPSILON_I_theta_l;
clear YPSILON_I_theta_l

% define in case there is no quadratic interpolation used
theta_ml=[];
YPSILON_I_theta_ml={};

if quad_interp==1
    theta_ml=theta_0-(theta_l-theta_0);
    YPSILON_I_theta_ml=cell(1,L);
    
    % create a waitbar to observe progress
    D = parallel.pool.DataQueue;
    h = waitbar(0,'Progress of parallel computations','Name','Please wait ...');
    p=1;
    afterEach(D, @nUpdateWaitbar);

    fprintf('\nCalculating interface modes at all minus support points in parallel (quadratic interpolation used)...\n');
    parfor k=1:L
        
        M_I_theta_ml=M_I_0;
        K_I_theta_ml=K_I_0;
        
        for l=1:n_theta
            M_I_theta_ml=M_I_theta_ml+M_I_j{l}*func_g{l}(theta_ml(l,k));
            K_I_theta_ml=K_I_theta_ml+K_I_j{l}*func_h{l}(theta_ml(l,k));
        end
        
        fprintf(['\nCalculating interface modes at support point -',num2str(k),' (quadratic interpolation used)...\n']);                
        [YPSILON_I_theta_ml{k},~]=interface_modes(K_I_theta_ml,M_I_theta_ml,dofdata,indata);
        
        send(D, k); % update waitbar
    end
    delete(h);
    clear M_I_theta_ml K_I_theta_ml 
end
 

submatdata_interp.theta_ml=theta_ml;
clear theta_ml
submatdata_interp.YPSILON_I_theta_ml=YPSILON_I_theta_ml;
clear YPSILON_I_theta_ml

% waitbar
function nUpdateWaitbar(~)
    waitbar(p/L, h);
    p = p + 1;
end

end

function [YPSILON,indata]=interface_modes(K_I,M_I,dofdata,indata)

n_IR=indata.n_IR;
method_eigf=indata.eigf.interface.method;
target_eigf=indata.eigf.interface.target;
max_eigf=indata.eigf.interface.max;
step_eigf=indata.eigf.interface.step;
init_eigf=indata.eigf.interface.init;

reduction_I=indata.reduction_I;

if reduction_I==5 % global with parametrization
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
    

    YPSILON=cell(1,N_I);
    for k=1:N_I
        fprintf(['    Calculating modes of interface ',num2str(k),'...\n']);
        [~,~,index]=intersect(I_dofs_l{k},I_dofs,'stable');
        M_Ill_l=M_I(index,index);
        K_Ill_l=K_I(index,index);
        [YPSILON{k},~,n_IR_l(k)]=keptmodes(K_Ill_l,M_Ill_l,method_eigf,n_IR_l(k),target_eigf(k),max_eigf(k),step_eigf(k),init_eigf(k));
        
        
        % Normalize YPSILON wrt the mass
        for l=1:size(YPSILON{k},2)
            YPSILON{k}(:,l)=YPSILON{k}(:,l)/sqrt(YPSILON{k}(:,l)'*M_Ill_l*YPSILON{k}(:,l));
        end
        
    end

    % update number of kept modes and dimensions of reduced matrices based on
    % the kept modes -> update input
    indata.n_IR_l=n_IR_l;
    n_IRL=sum(n_IR_l);
    indata.n_IRL=n_IRL;
    indata.n_DIL=indata.n_id+indata.n_IRL;
    
end

end