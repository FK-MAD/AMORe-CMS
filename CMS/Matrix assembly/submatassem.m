function [indata]=submatassem(dofdata,indata)

if ~isfile('mats_S.mat') % file does not exist
    fprintf('\n\n---Calculating sub-matrices---\n');
    
    [indata]=make_mats_S_k(dofdata,indata); % creates mats_S_k.mat (k=1,2,...)
    make_mats_S(indata); % creates mats_S.mat from mats_S_k.mat (k=1,2,...)
    make_mats(indata); % creates mats.mat from mats_S_k.mat (k=1,2,...) OR mats_S.mat (depending on indata.blkdiag_method)
else % file already exists
    fprintf('files mats_S_k.mat (k=1,2,...) already exist and will be used instead of creating new ones...\n\n');
    
    % find the number of kept normal modes used for each group
    N_S=indata.N_S;
    n_id_S_correct=zeros(1,N_S);
    for k=1:N_S
        filename=['mats_S_',num2str(k),'.mat'];
        mats_S_k=matfile(filename,'writable',false);
        n_id_S_correct(k)=size(mats_S_k.LAMBDA_id_S,2);
    end
    
    % number of kept fixed-interface normal modes (internal dofs) for each
    % group does not correspond to user's input || user wants to calculate
    % kept modes based on cutoff frequency -> n_id_S might become different 
    if ~isequal(indata.n_id_S,n_id_S_correct) || indata.eigf.group.method==1 
        warning(['The number of kept fixed-interface normal modes (internal dofs) for each group in mats_S.mat'...
        ' does not correspond to input vector n_id_S or method_eigf_S is set to 1 (kept modes based on cutoff frequency).']);
    
        fprintf('The input vector n_id_S is\n\n');
        disp(indata.n_id_S)
        
        fprintf('The vector n_id_S based on files mats_S_k.mat (k=1,2,...) is\n\n');
        disp(n_id_S_correct)    
        
        answer=input(['Enter "ok" to change input vector n_id_S to match that used to create files mats_S_k.mat (k=1,2,...).\n'...
            'Alternatively, press the Return key to update mats_S_k.mat (k=1,2,...) to correspond to input vector n_id_S or the specified cutoff frequencies.\n'...
            'To create new mats_S_k.mat (k=1,2,...), delete the existing files and re-run the program.\n']);
        if string(answer)=="ok"
            indata.n_id_S=n_id_S_correct;
            indata.n_id=sum(indata.n_id_S);
            indata.n_D=indata.n_id+dofdata.n_I;
            indata.n_DI=indata.n_id+indata.n_IR;
            indata.n_DIL=indata.n_id+indata.n_IRL;
            
%             mats_S=matfile('mats_S.mat','Writable',true); % create connection with file without loading in memory
%             update_mats_S(mats_S,indata);
%             mats=matfile('mats.mat','Writable',true); % create connection with file without loading in memory
%             update_mats(mats,indata);
        else
            [indata]=update_mats_S_k(dofdata,indata);
            
            mats_S=matfile('mats_S.mat','Writable',true); % create connection with file without loading in memory
            update_mats_S(mats_S,indata);
            mats=matfile('mats.mat','Writable',true); % create connection with file without loading in memory
            update_mats(mats,indata);
        end
    end
    
end

save('indata.mat','-struct','indata','-v7.3');

end