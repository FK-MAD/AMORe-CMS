function mats=make_mats(indata)

mats=struct; % define the variable where all matrices are stored as empty structure
save('mats.mat','-struct','mats','-v7.3'); % create the variable
mats=matfile('mats.mat','Writable',true); % create connection with variable without loading in memory

% make block diagonal matrices from mats_S_k.mat (k=1,2,...) and save them in file mats.mat 
mats.M_ii=blkdiag_DISC('M_ii_S',indata);
mats.M_ib=blkdiag_DISC('M_ib_S',indata);
mats.M_bb=blkdiag_DISC('M_bb_S',indata);

mats.K_ii=blkdiag_DISC('K_ii_S',indata);
mats.K_ib=blkdiag_DISC('K_ib_S',indata);
mats.K_bb=blkdiag_DISC('K_bb_S',indata);

mats.PSI_ib=blkdiag_DISC('PSI_ib_S',indata);
mats.PHI_id=blkdiag_DISC('PHI_id_S',indata);
mats.LAMBDA_id=blkdiag_DISC('LAMBDA_id_S',indata);

mats.M_ib_hat=blkdiag_DISC('M_ib_S_hat',indata);
mats.M_bb_hat=blkdiag_DISC('M_bb_S_hat',indata);
mats.K_bb_hat=blkdiag_DISC('K_bb_S_hat',indata);

mats.M_ib_tilde=blkdiag_DISC('M_ib_S_tilde',indata);
mats.F_bar=blkdiag_DISC('F_S_bar',indata);

mats.Properties.Writable = false; % to prevent further changes

end