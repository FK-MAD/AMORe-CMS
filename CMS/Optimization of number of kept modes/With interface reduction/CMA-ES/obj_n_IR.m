function Y=obj_n_IR(theta,mats,mats_I,dofdata,indata,LAMBDA)

    indata.eigf.interface.target=theta*4.6;

    
    fprintf('\n------\n');
    
    tic
    [matdata,indata]=optmatassem_global([],mats,mats_I,dofdata,indata);
    toc
    
    tic
    LAMBDA_reduced=eigs(matdata.K_D_reduced,matdata.M_D_reduced,length(LAMBDA),'smallestabs');
    toc    
    
    error=max(abs((sqrt(LAMBDA)'-sqrt(LAMBDA_reduced)')./sqrt(LAMBDA)')*100);
    
    disp(error)
    disp(theta)
    disp(indata.n_IR)
    
    Y=abs(error-1);%+norm(theta)%norm(indata.n_id_S)/norm(indata.n_id_S_store(1:indata.N_S));
    fprintf('------\n');
end