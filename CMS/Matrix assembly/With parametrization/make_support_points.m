function indata=make_support_points(indata,dofdata)

% define in case there is no parametrization with interface reduction 
indata.theta_l=[];
if indata.reduction_I==5 % global
    indata=make_support_points_global(indata);
elseif indata.reduction_I==6 % local
    indata=make_support_points_local(indata,dofdata); 
end
indata.L=size(indata.theta_l,2); % number of support points

end