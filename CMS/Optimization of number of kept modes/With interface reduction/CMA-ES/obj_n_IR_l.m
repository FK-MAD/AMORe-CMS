function Y=obj_n_IR_l(theta,mats,mats_I,dofdata,indata,LAMBDA)

    N_I=dofdata.N_I;

    indata.eigf.interface.target=theta*4.6*ones(1,N_I);

    
    fprintf('\n------\n');
    
    % to increase speed consider calculating K_II and M_II as full (if
    % memory allows for this) in optmatassem_local
    % K_II=full(YPSILON_Ill)'*full(K_I)*full(YPSILON_Ill);
    % M_II=full(YPSILON_Ill)'*full(M_I)*full(YPSILON_Ill);
    tic
    [matdata,indata]=optmatassem_local([],mats,mats_I,dofdata,indata);
    toc
    
    tic
    LAMBDA_reduced=eigs(matdata.K_D_reduced,matdata.M_D_reduced,length(LAMBDA),'smallestabs');
    toc    
    
    error=max(abs((sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)')*100);
    
    %disp(error)
    disp(abs((sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)')*100)
    disp(theta)
    disp(indata.n_IR_l)
    
    Y=abs(error-1);%+norm(theta)%norm(indata.n_id_S)/norm(indata.n_id_S_store(1:indata.N_S));
    fprintf('------\n');
end