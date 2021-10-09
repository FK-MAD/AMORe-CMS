function [group_S,S_0,S_j]=groupassem(model)

dim=model.geom('geom1').getSDim;
ngtag=mphtags(model.nodeGroup); % tags of node groups

if isempty(ngtag)
    error('Model must have at least one node group (e.g. theta0 containing all groups of sub-structures).');
end

% define in case all groups are independent of model parameters
S_j={};

for k=1:length(ngtag)
    nglabel=string(model.nodeGroup(ngtag{k}).label); % label of node group -> denotes which model parameter the contained groups depend on
    num_theta=str2double(erase(nglabel,"theta")); % the number of the dependent parameter extracted from the label
    seltag=string(model.nodeGroup(ngtag{k}).members); % tags of selections contained in node group
    num_group={}; % clear for every node group
    for l=1:length(seltag)
        sellabel=string(model.selection(seltag{l}).label); % label of selection of node group -> denotes the number of the group
        seldomain=model.selection(seltag{l}).entities(dim); % domains contained in selection
        num_group{l}=str2double(erase(sellabel,"group")); % number of all selections contained in node group
        group_S{num_group{l}}=seldomain'; % each cell has the domains contained in every group
    end  
    
    if num_theta==0
        S_0=[num_group{:}]; % domains independent of model parameters
    else
        S_j{num_theta}=[num_group{:}]; % each cell has the groups which depend on the respective model parameter
    end
        
    
end


end