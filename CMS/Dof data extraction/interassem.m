function I_dofs_l=interassem(dofdata,indata)

bound=dofdata.bound;
bound_dofs=dofdata.bound_dofs;
I_dofs=dofdata.I_dofs;
adj_I=dofdata.adj_I;
group_l=indata.group_l;

empty=[];
I_dofs_l={}; % define in case there are no boundaries (a single domain) 
remaining=I_dofs;
    
fprintf('\n\n');
k=0;
while ~isempty(remaining)
    k=k+1;
    fprintf(['Working on dofs of interface ',num2str(k),'...\n'])
    % find the indices of the boundaries contained in the selection of
    % adj_I specified by group_l. It is quaranteed that adjacent interfaces will be selected
    % together.
    [~,index_b]=intersect(bound,vertcat(adj_I{group_l{k}}));
    I_dofs_l{k}=intersect(remaining,vertcat(bound_dofs{index_b}),'sorted'); % dofs of each independent interface
    remaining=setdiff(remaining,I_dofs_l{k},'stable');
    if isempty(I_dofs_l{k})
        empty=[empty,k];
    end
end
    
I_dofs_l(empty)=[];

end