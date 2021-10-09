function Y=obj_n_id_S_squares(theta,mats,PSI_ib_S,M_ii_S,M_ib_S,K_ii_S,PHI_id_S_store,LAMBDA_id_S_store,tol,dofdata,indata,LAMBDA)

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
    
    error=abs(sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)';
    
    disp(max(error,[],'all'))
    disp(theta')
    disp(indata.n_id_S)
    
    M=0;
    if any(error>tol)
        M=1000;
    end
    
    Y=sum(theta.^2)+M;
    fprintf('------\n\n');
end
