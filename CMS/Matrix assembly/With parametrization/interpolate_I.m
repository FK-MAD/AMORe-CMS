function [mats_interp_k,submatdata_interp,indata]=interpolate_I(mats_interp_k,submatdata_interp,submatdata_M_I_K_I,dofdata,indata)
% mats_interp_k -> variable in RAM OR matfile. works with both
% submatdata_interp -> variable in RAM OR matfile. works with both
% submatdata_M_I_K_I -> variable in RAM OR matfile. works with both

theta_k=indata.theta_k;

% check if current sample point (theta_k) lies inside the convex hull of
% the support points. Different check for global and local reduction.
% Global needs all the parameters to lie in the convex hull.
% Local needs only the active parameter combinations to lie in the convex
% hull.
checkhull=check_inhull(indata);

% if not -> calculate the interface modes in the sample point and add the sample point to the support points
if checkhull==0 
    fprintf('\nCalculating interface modes in the sample point and updating support points...\n');
    submatdata_interp_new=update_interp(submatdata_M_I_K_I,dofdata,indata);
    
    submatdata_interp.YPSILON_I_theta_l(1,end+1)={submatdata_interp_new.YPSILON_I_theta_new};
    
    if ~isempty(submatdata_interp_new.theta_new_ml)
        submatdata_interp.YPSILON_I_theta_ml(1,end+1)={submatdata_interp_new.YPSILON_I_theta_new_ml};
        submatdata_interp.theta_ml(:,end+1)=submatdata_interp_new.theta_new_ml;
    end
    indata.theta_l(:,end+1)=theta_k;
    indata.L=size(indata.theta_l,2);
end

% interpolate. if the sample point was added to the support points, the
% interpolation yields the interface modes already calculated
mats_interp_k=matassem_interp(mats_interp_k,submatdata_interp,submatdata_M_I_K_I,dofdata,indata);  

end