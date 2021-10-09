function LAMBDA=errorviz(indata,matdata)

% number of modes to use in error calculation
num_modes=indata.num_modes;

% this loads the modes of the unreduced model. It calculates them if they
% do not exist.
if isfile('LAMBDA.mat')
    load('LAMBDA.mat','LAMBDA');
else
    warning('File "LAMBDA.mat" does not exist. Will calculate from unreduced model.')
    
    moddata=load('moddata.mat');
    [V,LAMBDA]=eigs(moddata.mats.Kc,moddata.mats.Ec,num_modes,'smallestabs');
    LAMBDA=diag(LAMBDA);
    for k=1:length(LAMBDA)
        V(:,k)=V(:,k)/sqrt(V(:,k)'*moddata.mats.Ec*V(:,k));
    end
    
    save('LAMBDA.mat','LAMBDA','-v7.3');
end

% keep all modes up to "num_modes"
LAMBDA=LAMBDA(1:num_modes);

[V_D,LAMBDA_D]=eigs(matdata.K_D_reduced,matdata.M_D_reduced,num_modes,'smallestabs');
LAMBDA_D=diag(LAMBDA_D);
for k=1:length(LAMBDA_D)
    V_D(:,k)=V_D(:,k)/sqrt(V_D(:,k)'*matdata.M_D_reduced*V_D(:,k));
end

error_D=abs(sqrt(LAMBDA)'-sqrt(LAMBDA_D)')./sqrt(LAMBDA)';
max_error=max(error_D);
max_mode=find(error_D==max_error);

figure
semilogy(error_D,'-bo');
xticks(1:20);
xlim([0 21]);
grid on
xlabel('Eigenfrequency Number','interpreter','latex');
ylabel('Fractional Error','interpreter','latex');
%text(max_mode,max_error,['Maximum error=',num2str(max_error)],'HorizontalAlignment','right');

end

