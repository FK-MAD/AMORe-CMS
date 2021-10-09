function [dofdata,indata]=dofassem(model,moddata,indata)

N_S=indata.N_S; % # of groups
group_S=indata.group_S; % domain indices for each group
null=moddata.mats.Null;
n_all=moddata.xmsh.ndofs;
dim=moddata.dim;

% find fixed dofs
fixed=find(sum(null,2)==0)-1; % fixed dofs correspond to all-zero rows of 'Null' matrix (-1 to match the 0-based numbering of dofs by COMSOL)

%% dof information for each sub-structure
dofdata_s=dofassem_s(moddata);

b_dofs_s=dofdata_s.b_dofs_s;
b_elem_s=dofdata_s.b_elem_s;
i_dofs_s=dofdata_s.i_dofs_s;
i_elem_s=dofdata_s.i_elem_s;
bound_s=dofdata_s.bound_s;
bound=dofdata_s.bound;
bound_dofs=dofdata_s.bound_dofs;
bound_elem=dofdata_s.bound_elem;

%% dof information for each group of sub-structures

fprintf('\n\n');
for k=1:N_S
    
    fprintf(['Working on dofs of component group ',num2str(k),'...\n'])
    
    all_bound_S=sort(vertcat(bound_s{group_S{k}})); % boundaries of all sub-structures in group
    % find the boundaries that appear exactly once; they become the boundaries of group
    inds=diff(all_bound_S)>0;  
    bound_S{k}=all_bound_S([true;inds]&[inds;true]); % boundaries of group 
    
    %[inds,~]=find(bound==bound_S{k}'); 
    [~,inds]=intersect(bound,bound_S{k}); % find the indices in 'bound' of the boundaries of group  
    b_dofs_S{k}=unique(vertcat(bound_dofs{inds}),'sorted'); % boundary dofs of group
    
    % find the dofs that appear between the sub-structures in group; they are added to the internall dofs    
    all_b_dofs_S=sort(vertcat(b_dofs_s{group_S{k}})); % boundary dofs of all sub-structures in group
    extra_i_dofs_S=unique(setdiff(all_b_dofs_S,b_dofs_S{k}),'sorted'); 
    
    i_dofs_S{k}=sort([vertcat(i_dofs_s{group_S{k}});extra_i_dofs_S]); % internall dofs of group

    % discard fixed dofs
    b_dofs_S{k}=setdiff(b_dofs_S{k},fixed);
    i_dofs_S{k}=setdiff(i_dofs_S{k},fixed);
    
end

dofdata.b_dofs_s=b_dofs_s;
dofdata.b_elem_s=b_elem_s;
dofdata.bound_s=bound_s;
dofdata.i_dofs_s=i_dofs_s;
dofdata.i_elem_s=i_elem_s;
dofdata.bound=bound;
dofdata.bound_elem=bound_elem;
dofdata.bound_dofs=bound_dofs;
dofdata.b_dofs_S=b_dofs_S;
dofdata.i_dofs_S=i_dofs_S;
dofdata.bound_S=bound_S;
bound_I=unique(vertcat(bound_S{:}));
dofdata.bound_I=bound_I;

%% pass adjacency information about boundaries

% for k=1:N_S
%     boundaries=bound_S{k};
%     m=1;
%     adj_I{m,k}=[];
%     for l=boundaries
%         adjacent=mphgetadj(model,'geom1','boundary','boundary',l);
%         common_1=intersect(boundaries,adjacent);
%         common_2=intersect(adj_I{m,k},common_1);
%         if ~isempty(common_1) && ~isempty(common_2) || isempty(adj_I{m,k})
%             adj_I{m,k}=unique(vertcat(adj_I{m,k},l,common_1));
%         elseif ~isempty(common_1) && isempty(common_2) && ~isempty(adj_I{m,k})
%             m=m+1;
%             adj_I{m,k}=unique(vertcat(l,common_1));
%         end
%     end
% end

remaining=bound_I;
l=0;
while ~isempty(remaining)
    l=l+1;
    adj_I{l}=[];
    adj_bound=remaining(1);
    while ~isempty(adj_bound)
        k=adj_bound(1);
        adj=mphgetadj(model,'geom1','boundary','boundary',k);
        adj_bound=unique(vertcat(adj_bound,setdiff(intersect(adj,bound_I),adj_I{l})));
        adj_I{l}=unique(vertcat(adj_I{l},adj_bound));
        adj_bound(adj_bound==k)=[];
    end
    remaining=setdiff(remaining,vertcat(adj_I{:}));
end

dofdata.adj_I=adj_I;
%%

b_dofs=vertcat(b_dofs_S{:}); % boundary dofs of all groups (containing duplicates)
i_dofs=vertcat(i_dofs_S{:}); % internal dofs of all groups
[I_dofs,~,tildeinds]=unique(b_dofs,'sorted'); % interface dofs of all independent interfaces (no duplicates)
n_dofs=[i_dofs;I_dofs]; % dofs of all groups

dofdata.b_dofs=b_dofs;
dofdata.i_dofs=i_dofs;
dofdata.I_dofs=I_dofs;
dofdata.n_dofs=n_dofs;

dofdata.n_b=length(b_dofs); % number of boundary dofs
dofdata.n_i=length(i_dofs); % number of internal dofs
dofdata.n_I=length(I_dofs); % number of interface dofs
dofdata.n=length(n_dofs);  % number of free dofs
dofdata.n_all=n_all; % number of free+fixed dofs

indata.n_D=indata.n_id+dofdata.n_I;

%% specify grouping of interfaces

accept=0;
indata.group_l=num2cell(1:length(adj_I));
while accept==0
    I_dofs_l=interassem(dofdata,indata);
    dofdata.I_dofs_l=I_dofs_l;
    dofdata.N_I=length(I_dofs_l); % number of interfaces
    
    % delete all interface selections to avoid possible overides
    seltag=mphtags(model.selection);
    for k=1:length(seltag)
        sellabel=string(model.selection(seltag{k}).label);
        if contains(sellabel,'interface')
            model.component("comp1").selection.remove(seltag{k});
        end
    end
    
    % create one explicit selection (in COMSOL) for each group of
    % interfaces. This is an alternative way to visualize interface
    % grouping and check if it is acceptable. 
    % NOTE: The visualization of the original selection is kept in the
    % first MATLAB figure of the loop. In COMSOL it is changed continually.
    % Refer to the first MATLAB figure to construct cell array "group_l".
    for k=1:length(I_dofs_l)
        seltag=['inter_sel',num2str(k)];
        model.component("comp1").selection.create(seltag, "Explicit");
        sel=model.component("comp1").selection(seltag);
        sellabel=['interface',num2str(k)];
        sel.label(sellabel);
        sel.geom(dim-1);
        sel.set(vertcat(adj_I{indata.group_l{k}}));
    end
    
    dofviz(moddata,dofdata,indata,'');
    
    correct=0;
    while correct==0
        answer=input(['\n\nPlease specify the grouping of interfaces. Enter a cell array to define "group_l".\n',...
            'Alternatively, if the current grouping is acceptable, press the Return key.\n\n']);
        if iscell(answer)
            if isempty(setdiff(1:length(dofdata.adj_I),[answer{:}])) && length([answer{:}])==length(dofdata.adj_I) % if answer contains all interfaces exactly once
                indata.group_l=answer;
                correct=1;
            else
                fprintf('\nCell array must contain all interfaces exactly once. Please try again.\n');
            end
        else
            accept=1;
            correct=1;
        end
    end
end

% I_dofs_l=interassem(dofdata,indata);
% dofdata.I_dofs_l=I_dofs_l;
% dofdata.N_I=length(I_dofs_l); % number of interfaces

%%

 %% test -> IT WORKS! IT ALLOWS LOCAL REDUCTION TO WORK WITHOUT RE-ORDERING OF Y
 I_dofs_test=vertcat(I_dofs_l{:});
[~,tildeinds_test]=ismember(b_dofs,I_dofs_test);
% [~,test_inds]=intersect(I_dofs_test,I_dofs);
 n_dofs_test=[i_dofs;I_dofs_test]; % dofs of all groups
 dofdata.n_dofs=n_dofs_test;
 dofdata.I_dofs=I_dofs_test;
 %%

T_tilde=speye(dofdata.n_I);
%dofdata.T_tilde=T_tilde(tildeinds,:); % transformation matrix such that b_dofs=T_tilde*I_dofs
%dofdata.T_tilde=T_tilde(tildeinds,test_inds);
 dofdata.T_tilde=T_tilde(tildeinds_test,:);

%[map,~]=find(setdiff(0:n_all-1,fixed)==n_dofs);
[~,~,map]=intersect(setdiff(0:n_all-1,fixed),n_dofs,'stable');
T_G=speye(dofdata.n);
dofdata.T_G=T_G(map,:); % transformation matrix such that (free dofs)=T_G*n_dofs
end