function out=blkdiag_DISC(var,indata)
method=indata.blkdiag_method;
N_S=indata.N_S;
if method==1 % use files mats_S_k.mat, read variables in a for-loop and build the block-diagonal incrementally
    out=sparse(0,0);
    for k=1:N_S
        filename=['mats_S_',num2str(k),'.mat'];
        %fprintf(['\nGetting matrix ',var,' of group #',num2str(k),' from ',filename,'...\n']);
        in=matfile(filename,'Writable',false); % create connection with file
        %fprintf('Adding diagonal block...\n');
        out=blkdiag(out,in.(var));
    end
elseif method==2 % use files mats_S_k.mat, read variables in a for-loop and build the block-diagonal once
    block=cell(1,N_S);
    for k=1:N_S
        filename=['mats_S_',num2str(k),'.mat'];
        %fprintf(['\nGetting matrix ',var,' of group #',num2str(k),' from ',filename,'...\n']);
        in=matfile(filename,'Writable',false); % create connection with file
        block{k}=in.(var);
    end
    %fprintf('Making block diagonal matrix...\n');
    out=blkdiag(block{1:end-1},sparse(block{end}));
else % use file mats_S.mat and build the block-diagonal once
    %fprintf(['\nGetting cell array of matrices ',var,' of all groups from mats_S.mat...\n']);
    in=matfile('mats_S.mat','Writable',false); % create connection with file
    block=in.(var);
    %fprintf('Making block diagonal matrix...\n');
    out=blkdiag(block{1:end-1},sparse(block{end}));
end
end