function varargout=diagassem(set,varargin)
    
for k=1:length(varargin{1})
    if ~ismember(k,set)
        for l=1:length(varargin)
            varargin{l}{k}=sparse(size(varargin{l}{k},1),size(varargin{l}{k},2));
        end
    end
end

varargout=cell(1,nargout);
for k=1:nargout
    varargout{k}=blkdiag(varargin{k}{1:end-1},sparse(varargin{k}{end})); % this way blkdiag returns sparse matrix
end
    
end