function moddata=modeldata(model,indata)
E=indata.E;
rho=indata.rho;
nu=indata.nu;
N_S=indata.N_S;

geom=model.geom('geom1');

moddata.dim=geom.getSDim; % dimension of model
moddata.N_s=geom.getNDomains; % # of sub-structures = # of domains that the geometry is partitioned in

if ~isempty(setdiff(1:moddata.N_s,[indata.group_S{:}])) || length([indata.group_S{:}])~=moddata.N_s % all domains are assigned || exactly once
    error('Check that all domains are assigned to groups and that each domain is assigned exactly once.');
end

moddata.adj=geom.getUpDown; % adjacent domains to every boundary (0 means there is no adjacent domain)

fprintf('\nGetting mesh information...\n');
[~,moddata.msh]=mphmeshstats(model); % mesh information

fprintf('\nGetting extended mesh information...\n');
moddata.xmsh=mphxmeshinfo(model); % extended mesh information

seltag=mphtags(model.selection);
for k=1:length(seltag)
    sellabel=string(model.selection(seltag{k}).label);
    num_theta(k)=str2double(erase(sellabel,"group")); % if NaN the group label is invalid
end

% delete all materials to avoid possible overides of selections
mattag=mphtags(model.material); % tags of existing materials
for k=1:length(mattag)    
    model.component("comp1").material.remove(mattag{k}); 
end


for k=1:N_S
    tag=['mat',num2str(k)];
%     if ~ismember(tag,string(mattag)) % if there is no material for sub-structure
%         model.component('comp1').material.create(tag); % create a material node
%     end
    model.component('comp1').material.create(tag); % create a material node
    mat=model.component('comp1').material(tag);
    % set the properties of sub-structure
    mat.propertyGroup('def').set('youngsmodulus',E(k));
    mat.propertyGroup('def').set('density',rho(k));
    mat.propertyGroup('def').set('poissonsratio',nu(k));
    [~,index_sel]=intersect(num_theta,k); % find the tag of the selection which is named 'groupk'
%     mat.selection.set([]); % clear the material of all selections
    mat.selection.named(seltag{index_sel}); % assign whole group (all domains of group) to the material
end

% get the eliminated stiffnes (Kc) and eliminated mass (Ec) matrices (free
% dofs only) along with the Null matrix (to find the fixed dofs)
% assembled using the solver node 'sol1'
%moddata.mats=mphmatrix(model,'sol1','out',{'K','E','Kc','Ec','Null'});
fprintf('\nGetting eliminated mass (Ec) and stiffness (Kc) matrices and matrix to find fixed dofs (Null)...\n');
moddata.mats=mphmatrix(model,'sol1','out',{'Kc','Ec','Null'}); 

if isempty(moddata.mats.Ec(moddata.mats.Ec~=0)) % check if mass matrix is empty
   error('Mass matrix is empty. Consider setting an eigenfrequency study.'); 
end

end