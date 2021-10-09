function mats_S=update_mats_S(mats_S,indata)
% mats_S -> variable in RAM OR matfile. works with both

N_S=indata.N_S;

% copy updated matrices from mats_S_k.mat(k=1,2,...) in cell arrays of mats_S.mat
fprintf('\nUpdating cell arrays of mats_S.mat...\n');
for k=1:N_S
    filename=['mats_S_',num2str(k),'.mat'];
    fprintf(['\nCopying data from file ',filename,' to mats_S.mat...\n']);
    mats_S_k=matfile(filename,'Writable',false); % create connection with variable without loading in memory
    
    mats_S.PHI_id_S(1,k)={mats_S_k.PHI_id_S};
    mats_S.LAMBDA_id_S(1,k)={mats_S_k.LAMBDA_id_S};
    
    mats_S.PHI_id_S_store(1,k)={mats_S_k.PHI_id_S_store};
    mats_S.LAMBDA_id_S_store(1,k)={mats_S_k.LAMBDA_id_S_store};

    mats_S.M_ib_S_hat(1,k)={mats_S_k.M_ib_S_hat};
    mats_S.F_S_bar(1,k)={mats_S_k.F_S_bar};
end

end