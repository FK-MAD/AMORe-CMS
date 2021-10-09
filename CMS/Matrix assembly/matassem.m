function [matdata,indata]=matassem(dofdata,indata)

reduction_I=indata.reduction_I;

mats=matfile('mats.mat','Writable',false); 
%mats=load('mats.mat'); 

% % if you want to save variables directly in matfiles -> slow but needs less
% % memory
% matdata=struct; % define the variable where all matrices are stored as empty structure
% save('matdata.mat','-struct','matdata','-v7.3'); % create the variable
% matdata=matfile('matdata.mat','Writable',true); % create connection with variable without loading in memory

% if you want to save variables in memory -> fast but needs more memory
 matdata=struct;

fprintf('\n\n---Assembling reduced matrices from sub-matrices---\n\n\n');
if reduction_I==1
    matdata=matassem_unreduced(matdata,mats,dofdata,indata);
elseif reduction_I==2
    [matdata,indata]=matassem_global(matdata,mats,dofdata,indata);
elseif reduction_I==3
    [matdata,indata]=matassem_local(matdata,mats,dofdata,indata);
else
%     % if you want to save variables directly in matfiles -> slow but needs less
%     % memory    
%     submatdata=struct; % define the variable where all matrices are stored as empty structure
%     save('submatdata.mat','-struct','submatdata','-v7.3'); % create the variable
%     submatdata=matfile('submatdata.mat','Writable',true); % create connection with variable without loading in memory    
    
    % if you want to save variables in memory -> fast but needs more memory
    submatdata=struct;
    
    mats_S=matfile('mats_S.mat','Writable',false); 
    %mats_S=load('mats_S.mat');
    if reduction_I==4       
        submatdata=submatassem_para_unreduced(submatdata,mats_S,mats,dofdata,indata);
        matdata=matassem_para_unreduced(matdata,submatdata,dofdata,indata);
    else
%         % if you want to save variables directly in matfiles -> slow but needs less
%         % memory
%         submatdata_interp=struct; % define the variable where all matrices are stored as empty structure
%         save('submatdata_interp.mat','-struct','submatdata_interp','-v7.3'); % create the variable
%         submatdata_interp=matfile('submatdata_interp.mat','Writable',true); % create connection with variable without loading in memory        
        
        % if you want to save variables in memory -> fast but needs more memory
        submatdata_interp=struct;
        
%         % if you want to save variables directly in matfiles -> slow but needs less
%         % memory
%         submatdata_M_I_K_I=struct; % define the variable where all matrices are stored as empty structure
%         save('submatdata_M_I_K_I.mat','-struct','submatdata_M_I_K_I','-v7.3'); % create the variable
%         submatdata_M_I_K_I=matfile('submatdata_M_I_K_I.mat','Writable',true); % create connection with variable without loading in memory      
        
        % if you want to save variables in memory -> fast but needs more memory
        submatdata_M_I_K_I=struct;
        
        submatdata_M_I_K_I=submatassem_M_I_K_I(submatdata_M_I_K_I,mats_S,dofdata,indata);
        [submatdata_interp,indata]=submatassem_interp(submatdata_interp,submatdata_M_I_K_I,dofdata,indata);
        
        theta_k=indata.theta_k;
        indata.theta_k=indata.theta_nom;
        mats_interp_nom=matassem_interp([],submatdata_interp,submatdata_M_I_K_I,dofdata,indata);
        
        if reduction_I==5
            submatdata=submatassem_para_global(submatdata,mats_S,mats,mats_interp_nom,dofdata,indata);
            
            indata.theta_k=theta_k;
            checkhull=check_inhull(indata);
            mats_interp_k=matassem_interp([],submatdata_interp,submatdata_M_I_K_I,dofdata,indata);
            matdata=matassem_para_global(matdata,submatdata,mats_interp_k,dofdata,indata);
        elseif reduction_I==6
            submatdata=submatassem_para_local(submatdata,mats_S,mats,mats_interp_nom,dofdata,indata);
            
            indata.theta_k=theta_k;
            checkhull=check_inhull(indata);
            mats_interp_k=matassem_interp([],submatdata_interp,submatdata_M_I_K_I,dofdata,indata);
            matdata=matassem_para_local(matdata,submatdata,mats_interp_k,dofdata,indata);
        end
        % in case you have used the option to save variables directly in matfile
        submatdata_interp.Properties.Writable = false; % to prevent further changes
        
        % save variables also in a matfile in case you have used the option
        % to save variables in memory
        % save('submatdata_interp.mat','-struct','submatdata_interp','-v7.3');
    end
    submatdata.Properties.Writable = false; % to prevent further changes
end

% in case you have used the option to save variables directly in matfile
matdata.Properties.Writable = false; % to prevent further changes

% save variables also in a matfile in case you have used the option
% to save variables in memory
% save('matdata.mat','-struct','matdata','-v7.3');

end