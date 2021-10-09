theta_0=[.5 .5 .5]';
theta_k=[.5 .5 1]';
theta_l=[1 0 0;0 1 0;0 0 1;.5 .4 1]';
keep=1;
theta_0=theta_0(1:keep,:);
theta_k=theta_k(1:keep,:);
theta_l=theta_l(1:keep,:);
L=size(theta_l,2);

[~,nearest]=min(vecnorm(theta_l-repmat(theta_k,1,size(theta_l,2)),2,1));

theta_q=theta_l(:,nearest);

xi_q_k=(theta_k-theta_0)'*(theta_q-theta_0)/norm(theta_q-theta_0)^2;

if isnan(xi_q_k)
    xi_q_k=1;
end

if nearest>L % if nearest point belongs to set theta(-l)
    xi_q_k=-xi_q_k; 
    nearest=nearest-L; % nearest must be between 1 and L when deleting the column of A that coresponds to theta_q
end

v_k=(theta_k-theta_0)-xi_q_k*(theta_q-theta_0);

A=theta_l-repmat(theta_0,1,L);
A(:,nearest)=[];

% [U,S,V]=svd(A,0);
% 
% tau_k=V*((U'*v_k)./diag(S));

tau_k=lsqminnorm(A,v_k);

xi_k=[tau_k(1:nearest-1);xi_q_k;tau_k(nearest:end)];