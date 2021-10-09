function [indata]=make_mats_S_k(dofdata,indata)

static=indata.static;

N_S=indata.N_S;

n_id_S=indata.n_id_S;

method_eigf=indata.eigf.group.method;
target_eigf=indata.eigf.group.target;
max_eigf=indata.eigf.group.max;
step_eigf=indata.eigf.group.step;
init_eigf=indata.eigf.group.init;

n_id_S_store=indata.n_id_S_store;

method_eigf_store=indata.eigf.group.method_store;
target_eigf_store=indata.eigf.group.target_store;
max_eigf_store=indata.eigf.group.max_store;
step_eigf_store=indata.eigf.group.step_store;
init_eigf_store=indata.eigf.group.init_store;

% create a waitbar to observe progress
D = parallel.pool.DataQueue;
h = waitbar(0,'Progress of parallel computations','Name','Please wait ...');
p=1;
afterEach(D, @nUpdateWaitbar);

fprintf('\nCalculating matrices for all groups in parallel and saving to files mats_S_k.mat (k=1,2,...)...\n');

parfor k=1:N_S
    filename=['mats_S_',num2str(k),'.mat'];
    mats_S_k=matfile(filename,'Writable',true); % create connection with variable without loading in memory

    %fprintf('        Saving...\n');   
    M_ii_S=mats_S_k.M_ii_S;
    M_ib_S=mats_S_k.M_ib_S;
    M_bb_S=mats_S_k.M_bb_S;
    K_ii_S=mats_S_k.K_ii_S;
    K_ib_S=mats_S_k.K_ib_S;
    K_bb_S=mats_S_k.K_bb_S;
    
    % interior partition of the interface constraint modes matrix
    %fprintf('    Calculating interface constrained modes Ψ...\n');
    PSI_ib_S=full(-K_ii_S\K_ib_S); % dense matrix (all elements non-zero) -> full (sparse makes it approximately two times larger)

    %fprintf('        Saving...\n');   
    mats_S_k.PSI_ib_S=PSI_ib_S;
    
    % matrices Φ_id and Λ_id containing the kept fixed-interface normal
    % modes
    %fprintf('    Calculating fixed-interface normal modes Φ...\n');
    % PHI_id_S contains kept modeshapes (all elements non-zero) -> full 
    [PHI_id_S,LAMBDA_id_S,n_id_S(k)]=keptmodes(K_ii_S,M_ii_S,method_eigf,n_id_S(k),target_eigf(k),max_eigf(k),step_eigf(k),init_eigf(k));
    LAMBDA_id_S=sparse(LAMBDA_id_S); % LAMBDA_id_S is diagonal and contains kept eigenvalues -> sparse 
    
    % Normalize PHI_id_S wrt the mass
    for l=1:n_id_S(k)
        PHI_id_S(:,l)=PHI_id_S(:,l)/sqrt(PHI_id_S(:,l)'*M_ii_S*PHI_id_S(:,l));
    end
    
    %fprintf('        Saving...\n');
    mats_S_k.PHI_id_S=PHI_id_S;
    mats_S_k.LAMBDA_id_S=LAMBDA_id_S;
    
    %% Store a big number of modes (more than n_id_S(k)) to use them without solving the eigenproblem when updating the matrices
    [PHI_id_S_store,LAMBDA_id_S_store,n_id_S_store(k)]=keptmodes(K_ii_S,M_ii_S,...
        method_eigf_store,n_id_S_store(k),target_eigf_store(k),max_eigf_store(k),step_eigf_store(k),init_eigf_store(k));
    LAMBDA_id_S_store=sparse(LAMBDA_id_S_store); % LAMBDA_id_S is diagonal and contains kept eigenvalues -> sparse 
    
    % Normalize PHI_id_S wrt the mass
    for l=1:n_id_S_store(k)
        PHI_id_S_store(:,l)=PHI_id_S_store(:,l)/sqrt(PHI_id_S_store(:,l)'*M_ii_S*PHI_id_S_store(:,l));
    end
    
    %fprintf('        Saving...\n');
    mats_S_k.PHI_id_S_store=PHI_id_S_store;
    mats_S_k.LAMBDA_id_S_store=LAMBDA_id_S_store;    
    
    %%
    
    % sub-matrices used in the assembly of reduced-order matrices (eq.
    % 1.25, 1.26, 1.27)
    %fprintf('    Calculating matrices in Eqs. (1.25), (1.26) and (1.27) and saving...\n');
    mats_S_k.M_ib_S_hat=PHI_id_S'*M_ii_S*PSI_ib_S+PHI_id_S'*M_ib_S; % dense (all elements non-zero) -> full
    mats_S_k.K_bb_S_hat=full(K_ib_S'*PSI_ib_S+K_bb_S); % dense (all elements non-zero) -> full
    
    K_ib_S=[];
    K_bb_S=[];

    mats_S_k.M_bb_S_hat=full((PSI_ib_S'*M_ii_S+M_ib_S')*PSI_ib_S+...
        PSI_ib_S'*M_ib_S+M_bb_S); % dense (all elements non-zero) -> full
 
    M_bb_S=[];
    
    % needed for the static correction
    % define in case there is no static correction
    mats_S_k.M_ib_S_tilde=[];
    mats_S_k.F_S_bar=[];
    if static==1
        % M_ib_S_tilde is dense (all elements non-zero) -> full
        %fprintf('    Calculating matrix in Eq. (1.34) for static correction and saving...\n');
        mats_S_k.M_ib_S_tilde=full(M_ib_S+M_ii_S*PSI_ib_S); % (1.34) using (1.9)
        M_ii_S=[];
        M_ib_S=[];
        PSI_ib_S=[];

        %fprintf('    Calculating the inverse of K_ii for static correction...\n');
        %K_ii_S_inv=K_ii_S\speye(size(K_ii_S));
        %F_S_bar=(K_ii_S^-1)-PHI_id_S*(LAMBDA_id_S^-1)*PHI_id_S'; % (1.36)     
        %F_S_bar=K_ii_S_inv-PHI_id_S*(LAMBDA_id_S^-1)*PHI_id_S'; % (1.36)  
        
        % F_S_bar is dense (all elements non-zero) -> full
        %fprintf('    Calculating matrix in Eq. (1.36) for static correction and saving...\n');   
        mats_S_k.F_S_bar=full(K_ii_S\(speye(size(K_ii_S))-K_ii_S*PHI_id_S*(LAMBDA_id_S^-1)*PHI_id_S')); % (1.36) 
        K_ii_S=[];
        PHI_id_S=[];
        LAMBDA_id_S=[];
    end
    send(D, k); % update waitbar
    mats_S_k.Properties.Writable = false; % to prevent further changes

end

delete(h);

% update number of kept modes and dimensions of reduced matrices based on
% the kept modes -> update input
indata.n_id_S=n_id_S;
n_id=sum(n_id_S);
indata.n_id=n_id;
indata.n_D=indata.n_id+dofdata.n_I;
indata.n_DI=indata.n_id+indata.n_IR;
indata.n_DIL=indata.n_id+indata.n_IRL;

indata.n_id_S_store=n_id_S_store;

% waitbar
function nUpdateWaitbar(~)
    waitbar(p/N_S, h);
    p = p + 1;
end

end