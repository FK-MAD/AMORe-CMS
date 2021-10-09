function mats=update_mats(mats,indata)
% mats -> variable in RAM OR matfile. works with both
fprintf('\nUpdating mats.mat...\n');

mats.PHI_id=blkdiag_DISC('PHI_id_S',indata);
mats.LAMBDA_id=blkdiag_DISC('LAMBDA_id_S',indata);
mats.M_ib_hat=blkdiag_DISC('M_ib_S_hat',indata);
mats.F_bar=blkdiag_DISC('F_S_bar',indata);

end