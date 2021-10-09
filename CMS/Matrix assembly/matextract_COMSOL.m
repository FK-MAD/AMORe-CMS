function matextract_COMSOL(model,dofdata,indata)

N_S=indata.N_S;

for k=1:N_S
    mats_S_k=struct; % define temporary variable (one for each group) where all matrices are stored as empty structure
    filename=['mats_S_',num2str(k),'.mat'];
    save(filename,'-struct','mats_S_k','-v7.3'); % create the variable
    mats_S_k=matfile(filename,'Writable',true); % create connection with variable without loading in memory
    
    fprintf(['\n\nExtracting mass and stiffness matrices of group ',num2str(k),' from COMSOL...\n']);
    [mats_S_k.M_ii_S,mats_S_k.M_ib_S,mats_S_k.M_bb_S,mats_S_k.K_ii_S,mats_S_k.K_ib_S,mats_S_k.K_bb_S]=...
        submatextract_COMSOL(model,dofdata,indata,k); % low density matrices -> all sparse  
    fprintf(['Saved partitions (internall and boundary) of matrices in ',filename,'\n']);
    mats_S_k.Properties.Writable = false; % to prevent further changes
end

end