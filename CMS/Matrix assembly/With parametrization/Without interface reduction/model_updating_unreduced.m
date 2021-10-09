function [matdata,indata]=model_updating_unreduced(dofdata,indata)

if ~isfile('submatdata_4.mat') 
    submatdata=matfile('submatdata_4.mat','Writable',true);
    
    mats=matfile('mats.mat','Writable',false);
    mats_S=matfile('mats_S.mat','Writable',false);
    
    submatdata=submatassem_para_unreduced(submatdata,mats_S,mats,dofdata,indata);
    submatdata.Properties.Writable = false; % to prevent further changes
end
submatdata=load('submatdata_4.mat');

matdata=matassem_para_unreduced([],submatdata,dofdata,indata);
    
% l=0;
% for k=1:10 % dummy loop to test speed
%     l=l+1;
%     fprintf(['\nIteration #',num2str(l),'\n']);
%     indata.theta_k=normrnd(1,.05,4,1); % change the sample point      
%     
%     tic
%     matdata=matassem_para_unreduced([],submatdata,dofdata,indata);
%     toc
% end

end