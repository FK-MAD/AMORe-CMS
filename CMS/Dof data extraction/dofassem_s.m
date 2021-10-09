function dofdata_s=dofassem_s(moddata)

dim=moddata.dim;
N_s=moddata.N_s; % # of sub-structures
adj=moddata.adj;
msh=moddata.msh;
xmsh=moddata.xmsh;

fprintf('\n\n');

if dim==2
    %% mesh information
    % find the column where information about 'edge' elements are stored -> 1D elements used on the boundaries
    %[~,edgindex]=find(string(msh.types)=='edg');
    [~,edgindex]=intersect(string(msh.types),"edg",'stable');
    
    % find the 2D elements used
%     types=["tri","quad"]';
%     [inds,mshtypeindex_2D]=find(string(msh.types)==types);
%     mshtype_2D=types(inds);

types=["tri","quad"];
[mshtype_2D,mshtypeindex_2D]=intersect(string(msh.types),types,'stable');
    
    %% dof information for each sub-structure
    for k=1:N_s
        fprintf(['Working on dofs of geometrical domain ',num2str(k),'...\n']);
        [~,bound_s{k}]=find(adj(1,:)==k & adj(2,:)~=0 |...
            adj(2,:)==k & adj(1,:)~=0); % find the common boundaries of the domain k
        
        [b_elem_s{k},~]=find(msh.elementity{edgindex}==bound_s{k}); % find the edge elements belonging to all common boundaries
        bound_s{k}=bound_s{k}';
        b_dofs_s{k}=xmsh.elements.edg.dofs(:,b_elem_s{k}); % find the dofs belonging to edge elements of all common boundaries
        b_dofs_s{k}=unique(b_dofs_s{k}(:),'sorted'); % discard duplicate dofs and sort -> boundary dofs in vector u_b of sub-structure k
        
        i_dofs_s{k}=[];
        for l=1:length(mshtype_2D) % for every element type
            i_elem_s{l,k}=find(msh.elementity{mshtypeindex_2D(l)}==k); % find the elements belonging to domain k
            i_dofs_s_elem=xmsh.elements.(mshtype_2D(l)).dofs(:,i_elem_s{l,k}); % find the dofs belonging to elements of domain k
            i_dofs_s{k}=[i_dofs_s{k};i_dofs_s_elem(:)];
        end
        i_dofs_s{k}=unique(setdiff(i_dofs_s{k},b_dofs_s{k}),'sorted'); % discard boundary dofs and duplicate dofs and sort -> internal dofs in vector u_i of sub-structure k
    end
    
    bound=unique(vertcat(bound_s{:}),'sorted');
    
    % define in case there are no boundaries (a single domain)
    bound_elem={};
    bound_dofs={};
    fprintf('\n\n');
    for k=1:length(bound)
        fprintf(['Working on dofs of boundary ',num2str(k),'...\n']);
        [bound_elem{k},~]=find(msh.elementity{edgindex}==bound(k)); % find the edge elements belonging to boundary k
        bound_dofs{k}=xmsh.elements.edg.dofs(:,bound_elem{k}); % find the dofs belonging to edge elements of boundary k
        bound_dofs{k}=unique(bound_dofs{k}(:),'sorted'); % discard duplicate dofs and sort -> dofs of boundary k
    end
else
    % find the 2D elements used
    %types=["tri","quad"]';
    %[inds,mshtypeindex_2D]=find(string(msh.types)==types);
    types=["tri","quad"];
    [mshtype_2D,mshtypeindex_2D]=intersect(string(msh.types),types,'stable');
    
    %mshtype_2D=types(inds);
    
    % find the 3D elements used
%     types=["tet","pyr","prism","hex"]';
%     [inds,mshtypeindex_3D]=find(string(msh.types)==types);
%     mshtype_3D=types(inds);

types=["tet","pyr","prism","hex"];
[mshtype_3D,mshtypeindex_3D]=intersect(string(msh.types),types,'stable');

    
    %% dof information for each sub-structure
    for k=1:N_s
        fprintf(['Working on dofs of geometrical domain ',num2str(k),'...\n']);
        [~,bound_s{k}]=find(adj(1,:)==k & adj(2,:)~=0 |...
            adj(2,:)==k & adj(1,:)~=0); % find the common boundaries of the domain k
        
        b_dofs_s{k}=[];
        for l=1:length(mshtype_2D) % for every 2D element type
            [b_elem_s{l,k},~]=find(msh.elementity{mshtypeindex_2D(l)}==bound_s{k}); % find the elements belonging to all common boundaries
            b_dofs_s_elem=xmsh.elements.(mshtype_2D(l)).dofs(:,b_elem_s{l,k}); % find the dofs belonging to elements of all common boundaries
            b_dofs_s{k}=[b_dofs_s{k};b_dofs_s_elem(:)];
        end
        bound_s{k}=bound_s{k}';
        b_dofs_s{k}=unique(b_dofs_s{k}(:),'sorted'); % discard duplicate dofs and sort-> boundary dofs in vector u_b of sub-structure k
        
        i_dofs_s{k}=[];
        for l=1:length(mshtype_3D) % for every 3D element type
            i_elem_s{l,k}=find(msh.elementity{mshtypeindex_3D(l)}==k); % find the elements belonging to domain k
            i_dofs_s_elem=xmsh.elements.(mshtype_3D(l)).dofs(:,i_elem_s{l,k}); % find the dofs belonging to elements of domain k
            i_dofs_s{k}=[i_dofs_s{k};i_dofs_s_elem(:)];
        end
        i_dofs_s{k}=unique(setdiff(i_dofs_s{k},b_dofs_s{k}),'sorted'); % discard boundary dofs and duplicate dofs and sort -> internal dofs in vector u_i of sub-structure k
    end
    
    bound=unique(vertcat(bound_s{:}),'sorted');
    
    % define in case there are no boundaries (a single domain)
    bound_elem={};
    bound_dofs={};
    fprintf('\n\n');
    for k=1:length(bound)
        fprintf(['Working on dofs of boundary ',num2str(k),'...\n']);
        bound_dofs{k}=[];
        for l=1:length(mshtype_2D) % for every 2D element type
            [bound_elem{l,k},~]=find(msh.elementity{mshtypeindex_2D(l)}==bound(k)); % find the elements belonging to boundary k
            bound_dofs_elem=xmsh.elements.(mshtype_2D(l)).dofs(:,bound_elem{l,k}); % find the dofs belonging to elements of boundary k
            bound_dofs{k}=[bound_dofs{k};bound_dofs_elem(:)];
        end
        bound_dofs{k}=unique(bound_dofs{k}(:),'sorted'); % discard duplicate dofs and sort -> dofs of boundary k
    end
    
end

dofdata_s.b_dofs_s=b_dofs_s;
dofdata_s.b_elem_s=b_elem_s;
dofdata_s.i_dofs_s=i_dofs_s;
dofdata_s.i_elem_s=i_elem_s;
dofdata_s.bound_s=bound_s;
dofdata_s.bound=bound;
dofdata_s.bound_dofs=bound_dofs;
dofdata_s.bound_elem=bound_elem;

end