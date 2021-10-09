function [mats,PHI_id_S_store_sliced,LAMBDA_id_S_store_sliced,indata]=...
    update_mats_RAM(mats,PSI_ib_S_sliced,M_ii_S_sliced,M_ib_S_sliced,...
    K_ii_S_sliced,PHI_id_S_store_sliced,LAMBDA_id_S_store_sliced,dofdata,indata)
% input & output variables -> all in RAM

static=indata.static;

N_S=indata.N_S;

n_id_S=indata.n_id_S;
method_eigf=indata.eigf.group.method;
target_eigf=indata.eigf.group.target;
max_eigf=indata.eigf.group.max;
step_eigf=indata.eigf.group.step;
init_eigf=indata.eigf.group.init;

parfor k=1:N_S
    %fprintf(['\n\nUpdating matrices of group ',num2str(k),'...\n']);
    PSI_ib_S=PSI_ib_S_sliced{k};
    M_ii_S=M_ii_S_sliced{k};
    M_ib_S=M_ib_S_sliced{k};
    K_ii_S=K_ii_S_sliced{k};
    
    PHI_id_S=PHI_id_S_store_sliced{k};
    LAMBDA_id_S=LAMBDA_id_S_store_sliced{k};
    
    % Get the needed modes from the stored ones. If they are not enough
    % solve the eigenproblem
    if method_eigf==0 && length(LAMBDA_id_S)>=n_id_S(k) % contain at least as many modes as requested
        PHI_id_S=PHI_id_S(:,1:n_id_S(k));
        LAMBDA_id_S=LAMBDA_id_S(1:n_id_S(k),1:n_id_S(k));
    elseif method_eigf==1 && max(LAMBDA_id_S,[],'all')>=(target_eigf(k)*2*pi)^2 % contain at least the modes corresponding to requested cutoff frequency
        kept=sum(diag(LAMBDA_id_S)<(target_eigf(k)*2*pi)^2)+1;
        n_id_S(k)=kept;
        PHI_id_S=PHI_id_S(:,1:kept);
        LAMBDA_id_S=LAMBDA_id_S(1:kept,1:kept);
    else % stored modes are not enough -> must calculate new eigenproblem
        warning(['Not enough modes stored for group ',num2str(k),'. Solving the eigenproblem...']);
        % matrices Φ_id and Λ_id containing the kept fixed-interface normal
        % modes -> need updating
        % fprintf('    Updating fixed-interface normal modes Φ...\n');
        % PHI_id_S contains kept modeshapes (all elements non-zero) -> full 
        [PHI_id_S,LAMBDA_id_S,n_id_S(k)]=keptmodes(K_ii_S,M_ii_S,method_eigf,n_id_S(k),target_eigf(k),max_eigf(k),step_eigf(k),init_eigf(k));
        LAMBDA_id_S=sparse(LAMBDA_id_S); % LAMBDA_id_S is diagonal and contains kept eigenvalues -> sparse 

        % Normalize PHI_id_S wrt the mass
        for l=1:n_id_S(k)
            PHI_id_S(:,l)=PHI_id_S(:,l)/sqrt(PHI_id_S(:,l)'*M_ii_S*PHI_id_S(:,l));
        end
        
        % update storage since the calculated modes are more than those
        % stored
        PHI_id_S_store_sliced{k}=PHI_id_S;
        LAMBDA_id_S_store_sliced{k}=LAMBDA_id_S;
    end
    
    %fprintf('        Saving...\n');
    PHI_id_S_sliced{k}=PHI_id_S;
    LAMBDA_id_S_sliced{k}=LAMBDA_id_S;
    
    % sub-matrix in Eq. (1.25) needs updating
    %fprintf('    Updating matrix in Eq. (1.25) and saving...\n');
    M_ib_S_hat_sliced{k}=PHI_id_S'*M_ii_S*PSI_ib_S+PHI_id_S'*M_ib_S; % dense (all elements non-zero) -> full
    PSI_ib_S=[];
    M_ii_S=[];
    M_ib_S=[];
    
    % needed for the static correction
    F_S_bar_sliced{k}=[];
    if static==1
        % F_S_bar needs updating
        % F_S_bar is dense (all elements non-zero) -> full
        %fprintf('    Calculating matrix in Eq. (1.36) for static correction and saving...\n');   
%         filename_dest=['mats_S_',num2str(k),'.mat']; % write to file (already created) and not into memory due to large size of variable
%         mats_S_dest=matfile(filename_dest,'Writable',true); % create connection with destination variable without loading in memory
%         mats_S_dest.F_S_bar=full(K_ii_S\(speye(size(K_ii_S))-K_ii_S*PHI_id_S*(LAMBDA_id_S^-1)*PHI_id_S')); % (1.36) 
        
        F_S_bar_sliced{k}=full(K_ii_S\(speye(size(K_ii_S))-K_ii_S*PHI_id_S*(LAMBDA_id_S^-1)*PHI_id_S')); % (1.36)
        K_ii_S=[];
        PHI_id_S=[];
        LAMBDA_id_S=[]; 
    end
end

% make block diagonals from temporary variables and pass them in mats
mats.PHI_id=blkdiag(PHI_id_S_sliced{1:end-1},sparse(PHI_id_S_sliced{end}));
clear PHI_id_S_sliced
mats.LAMBDA_id=blkdiag(LAMBDA_id_S_sliced{1:end-1},sparse(LAMBDA_id_S_sliced{end}));
clear LAMBDA_id_S_sliced
mats.M_ib_hat=blkdiag(M_ib_S_hat_sliced{1:end-1},sparse(M_ib_S_hat_sliced{end}));
clear M_ib_S_hat_sliced
mats.F_bar=blkdiag(F_S_bar_sliced{1:end-1},sparse(F_S_bar_sliced{end}));
clear F_S_bar_sliced

% update number of kept modes and dimensions of reduced matrices based on
% the kept modes -> update input
indata.n_id_S=n_id_S;
n_id=sum(n_id_S);
indata.n_id=n_id;
indata.n_D=indata.n_id+dofdata.n_I;
indata.n_DI=indata.n_id+indata.n_IR;
indata.n_DIL=indata.n_id+indata.n_IRL;

end