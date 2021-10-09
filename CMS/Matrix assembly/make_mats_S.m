function mats_S=make_mats_S(indata)

N_S=indata.N_S;

mats_S=struct; % define the variable where all matrices are stored as empty structure
save('mats_S.mat','-struct','mats_S','-v7.3'); % create the variable
mats_S=matfile('mats_S.mat','Writable',true); % create connection with variable without loading in memory

% copy matrices from mats_S_k.mat(k=1,2,...) in cell arrays of mats_S.mat
for k=1:N_S
    filename=['mats_S_',num2str(k),'.mat'];
    fprintf(['\nCopying data from file ',filename,' to mats_S.mat...\n']);
    mats_S_k=matfile(filename,'Writable',false); % create connection with variable without loading in memory
    
    mats_S.M_ii_S(1,k)={mats_S_k.M_ii_S};
    mats_S.M_ib_S(1,k)={mats_S_k.M_ib_S};
    mats_S.M_bb_S(1,k)={mats_S_k.M_bb_S};
    mats_S.K_ii_S(1,k)={mats_S_k.K_ii_S};
    mats_S.K_ib_S(1,k)={mats_S_k.K_ib_S};
    mats_S.K_bb_S(1,k)={mats_S_k.K_bb_S};
    
    mats_S.PSI_ib_S(1,k)={mats_S_k.PSI_ib_S};
    
    mats_S.PHI_id_S(1,k)={mats_S_k.PHI_id_S};
    mats_S.LAMBDA_id_S(1,k)={mats_S_k.LAMBDA_id_S};
    
    mats_S.PHI_id_S_store(1,k)={mats_S_k.PHI_id_S_store};
    mats_S.LAMBDA_id_S_store(1,k)={mats_S_k.LAMBDA_id_S_store};
    
    mats_S.M_ib_S_hat(1,k)={mats_S_k.M_ib_S_hat};
    mats_S.M_bb_S_hat(1,k)={mats_S_k.M_bb_S_hat};
    mats_S.K_bb_S_hat(1,k)={mats_S_k.K_bb_S_hat};

    mats_S.M_ib_S_tilde(1,k)={mats_S_k.M_ib_S_tilde};
    mats_S.F_S_bar(1,k)={mats_S_k.F_S_bar};
    
    %delete(filename); % delete mats_S_k.mat
end

mats_S.Properties.Writable = false; % to prevent further changes

end