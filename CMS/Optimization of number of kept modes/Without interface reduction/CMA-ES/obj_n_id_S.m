function Y=obj_n_id_S(theta,mats,PSI_ib_S,M_ii_S,M_ib_S,K_ii_S,PHI_id_S_store,LAMBDA_id_S_store,tol,dofdata,indata,LAMBDA)

    indata.eigf.group.target=theta'*2000;
    
    fprintf('\n------\n');
   
%     tic
%     [mats,indata]=update_mats_HYBRID(mats,dofdata,indata);
%     toc
    
    tic
    [mats,~,~,indata]=...
    update_mats_RAM(mats,PSI_ib_S,M_ii_S,M_ib_S,...
    K_ii_S,PHI_id_S_store,LAMBDA_id_S_store,dofdata,indata);
    toc
    
    tic
    matdata=matassem_unreduced([],mats,dofdata,indata);
    toc
    
    tic
    LAMBDA_reduced=eigs(matdata.K_D_reduced,matdata.M_D_reduced,length(LAMBDA),'smallestabs');
    toc    
    
    error=max(abs(sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)');
    
    disp(error)
    disp(theta')
    disp(indata.n_id_S)
    
    Y=(error-tol)^2+sum(theta.^2)*10^-13;%+norm(theta)%norm(indata.n_id_S)/norm(indata.n_id_S_store(1:indata.N_S));
    fprintf('------\n');
end
