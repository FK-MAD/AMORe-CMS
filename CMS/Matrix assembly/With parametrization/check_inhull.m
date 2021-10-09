function checkhull=check_inhull(indata)

theta_k=indata.theta_k;
theta_l=indata.theta_l;
reduction_I=indata.reduction_I;

% check if current sample point (theta_k) lies inside the convex hull of
% the support points
if reduction_I==5 % global
    % check all parameters at once
    checkhull=inhull(theta_k',theta_l',[],1.e-13*mean(abs(theta_l),'all'));
    if checkhull==0
        warning('Sample point does not lie in the convex hull of the support points.');
    end 
elseif reduction_I==6 % local
    
    kept_param_unique=indata.kept_param_unique;
    
    % check only the active parameter combinations one-by-one
    k=1;
    out_of_hull=false;
    while k<=length(kept_param_unique) && out_of_hull==false
        params=kept_param_unique{k};
        if length(params)>=2 % for at least 2 dimensions
            checkhull=inhull(theta_k(params)',theta_l(params,:)',[],1.e-13*mean(abs(theta_l(params,:)),'all'));
        else % for 1 dimension
            if min(theta_l(params,:))<=theta_k(params) && theta_k(params)<=max(theta_l(params,:))
                checkhull=true;
            else
                checkhull=false;
            end
        end
        if checkhull==0
            warning(['Sample point for parameter combination [',num2str(params),'] does not lie in the convex hull of the support points.']);
            out_of_hull=true;
        end
        k=k+1;
    end
end


end