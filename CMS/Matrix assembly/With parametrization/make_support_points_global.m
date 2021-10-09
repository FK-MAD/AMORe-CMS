function indata=make_support_points_global(indata)

scatter_theta_l=indata.scatter_theta_l;
method_theta_l=indata.method_theta_l;
theta_0=indata.theta_0;
n_theta=indata.n_theta;

if method_theta_l=="simplex"
    indata.theta_l=scatter_theta_l.*simplex_coordinates(n_theta)+theta_0;
elseif method_theta_l=="hypercube"   
    % find all possible sign combinations for the given dimension
    sign_pairs=repmat({[1 -1]},1,n_theta);
    sign_combinations=allcomb(sign_pairs{:});
    
    % fill the rows corresponding to the parameters with the scaled vertex coordinates
    indata.theta_l=scatter_theta_l.*sign_combinations'.*.5.*ones(n_theta,size(sign_combinations,1))+theta_0;
end



end