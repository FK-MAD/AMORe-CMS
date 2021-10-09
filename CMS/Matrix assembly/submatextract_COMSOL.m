function [M_ii_S,M_ib_S,M_bb_S,K_ii_S,K_ib_S,K_bb_S]=submatextract_COMSOL(model,dofdata,indata,k)
E=indata.E;
rho=indata.rho;
nu=indata.nu;
N_S=indata.N_S;
i_dofs_S=dofdata.i_dofs_S;
b_dofs_S=dofdata.b_dofs_S;

mat=model.component('comp1').material(['mat',num2str(k)]);
mat.propertyGroup('def').set('youngsmodulus',E(k));
mat.propertyGroup('def').set('density',rho(k));
mat.propertyGroup('def').set('poissonsratio',nu(k));
        
fprintf('    Setting E and ρ of all other groups to zero...\n');
rest=[1:k-1 k+1:N_S]; % all groups except k -> the properties of their domains are set to zero
for l=rest
    mat=model.component('comp1').material(['mat',num2str(l)]);
    mat.propertyGroup('def').set('youngsmodulus',0);
    mat.propertyGroup('def').set('density',0);
end

fprintf('    Getting mass (E) and stiffness (K) matrices...\n');
mats=mphmatrix(model,'sol1','out',{'K','E'});

% assembly of mass matrix and its partitions
M_ii_S=mats.E(i_dofs_S{k}+1,i_dofs_S{k}+1); % +1 because numbering of dofs is 0-based by COMSOL
M_ib_S=mats.E(i_dofs_S{k}+1,b_dofs_S{k}+1);
M_bb_S=mats.E(b_dofs_S{k}+1,b_dofs_S{k}+1);

% assembly of stiffness matrix and its partitions
K_ii_S=mats.K(i_dofs_S{k}+1,i_dofs_S{k}+1);
K_ib_S=mats.K(i_dofs_S{k}+1,b_dofs_S{k}+1);
K_bb_S=mats.K(b_dofs_S{k}+1,b_dofs_S{k}+1);


end