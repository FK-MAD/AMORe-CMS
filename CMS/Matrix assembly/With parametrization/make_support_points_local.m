function indata=make_support_points_local(indata,dofdata)

scatter_theta_l=indata.scatter_theta_l;
method_theta_l=indata.method_theta_l;
theta_0=indata.theta_0;
n_theta=indata.n_theta;
N_I=dofdata.N_I;

% n_theta=6;
% theta_0=0*ones(n_theta,1);
% scatter_theta_l=.1*ones(n_theta,1);
% scatter_theta_l(3)=.3;
% scatter_theta_l(4)=.4;
% 
% indata.n_theta=n_theta;
% indata.theta_0=theta_0;
% indata.scatter_theta_l=scatter_theta_l;


% find the parameters that each interface depends on
% e.g. for a model with 4 interfaces and 4 parameters
% interface 1 -> parameters 1 and 2
% interface 2 -> parameters 1,2 and 4
% interface 3 -> parameter 3
% interface 4 -> parameters 1 and 2
% kept_param={[1,2],[1,2,4],[3],[1,2]};

%kept_param={[1,2],[2,3],[3,4],[4,5],[5,6],[6,7],[7,8],[8,9],[9,10]}; % ONLY TEMPORARY! COMMENT IT!!
kept_param=cell(1,N_I); % ONLY TEMPORARY! UNCOMMENT IT!!
dim=zeros(1,N_I);
for k=1:N_I
    kept_param{k}=find_kept_param(dofdata,indata,k); %TEMPORARY FOR LOCAL VS GLOBAL!!! UNCOMMENT IT!!
    dim(k)=length(kept_param{k});
end

% kept_param={[1,2,6],[1,4,6],[4,5,6],[1,4],[1,3],2,6,1,3,4};
% dim=[3,3,3,2,2,1,1,1,1,1];

% reorder based on descending dimensions
% e.g. kept_param={[1,2,4],[1,2],[1,2],[3]}
[~,ordering]=sort(dim,'descend');
kept_param=kept_param(ordering);
dim=dim(ordering);

% keep only the unique parameter combination in each dimension
% e.g. kept_param={[1,2,4],[1,2],[3]}
[dim_unique,dim_index]=unique(dim,'stable');
kept_param_unique={};
l=1;
for k=1:length(dim_index)
    if k<length(dim_index)
        same_dim_positions=dim_index(k):dim_index(k+1)-1;
    else
        same_dim_positions=dim_index(k):length(kept_param);
    end
    same_dim_param=vertcat(kept_param{same_dim_positions});
    [~,unique_same_dim_param]=unique(same_dim_param,'rows','stable');
    kept_param_unique=[kept_param_unique kept_param(unique_same_dim_param+dim_index(k)-1)];
    unique_same_dim_positions{k}=l:l+length(unique_same_dim_param)-1;
    l=l+length(unique_same_dim_positions{k});
end

indata.kept_param_unique=kept_param_unique;
n_unique=length(kept_param_unique);

for k=1:length(unique_same_dim_positions) % loop for each dimension e.g. 3D, 2D, 1D
    best_combos{k}=find_best_combo(kept_param_unique,unique_same_dim_positions{k},dim_unique(k));
    inputs{k}=1:size(best_combos{k},3);
end

all_combos=allcomb(inputs{:});
        
if method_theta_l=="simplex"        
        
    theta_l_simplex=cell(1,n_unique);

        % make the simplex for each unique parameter combination
        for k=1:n_unique
            param_k=kept_param_unique{k};
            dim_k=length(param_k);
            % make a cell for each unique combination of parameters
            theta_l_simplex{k}=zeros(n_theta,dim_k+1);
            % fill the rows corresponding to the parameters with the scaled vertex coordinates
            theta_l_simplex{k}(param_k,:)=scatter_theta_l(param_k).*simplex_coordinates(dim_k)+theta_0(param_k);
        end
       
        
        theta_l_dim=cell(1,size(all_combos,1));
        theta_l_merged=cell(1,size(all_combos,1));
        num_points=zeros(1,size(all_combos,1));
        combo_cell=cell(1,length(best_combos));
        
        for k=1:size(all_combos,1)
            for l=1:length(best_combos)
                combo_cell{l}=best_combos{l}(:,:,all_combos(k,l));
            end
            
            theta_l_dim{k}=make_theta_l_dim(theta_l_simplex,indata,kept_param_unique,unique_same_dim_positions,combo_cell);
            theta_l_merged{k}=make_theta_l_merged(theta_l_dim{k},n_theta,unique_same_dim_positions);
            num_points(k)=size(theta_l_merged{k},2);
            %disp(num_points(k))
        end
        
        [~,best_index]=min(num_points);
        theta_l=theta_l_merged{best_index};
           
    
elseif method_theta_l=="hypercube"
    
    theta_l_hypercube=cell(1,n_unique);

    % make the hypercube for each unique parameter combination
    for k=1:n_unique
        param_k=kept_param_unique{k};
        dim_k=length(param_k);
        % make a cell for each unique combination of parameters
        theta_l_hypercube{k}=zeros(n_theta,2^dim_k);
        
        % find all possible sign combinations for the given dimension
        sign_pairs=repmat({[1 -1]},1,dim_k);
        sign_combinations=allcomb(sign_pairs{:});

        % fill the rows corresponding to the parameters with the scaled vertex coordinates
        theta_l_hypercube{k}(param_k,:)=scatter_theta_l(param_k).*sign_combinations'.*.5.*ones(dim_k,size(sign_combinations,1))+theta_0(param_k);
    end
    
    
    
    for k=1:size(all_combos,1)
        for l=1:length(best_combos)
            combo_cell{l}=best_combos{l}(:,:,all_combos(k,l));
        end
        
        theta_l_dim_hypercube{k}=make_theta_l_dim(theta_l_hypercube,indata,kept_param_unique,unique_same_dim_positions,combo_cell);
        theta_l_merged_hypercube{k}=make_theta_l_merged(theta_l_dim_hypercube{k},n_theta,unique_same_dim_positions);
        num_points_hypercube(k)=size(theta_l_merged_hypercube{k},2);
    end
    
    [~,best_index_hypercube]=min(num_points_hypercube);
    theta_l=theta_l_merged_hypercube{best_index_hypercube};

end

% fill zero elements of each row with values between the minimum and
% maximum element of the row
for k=1:n_theta
    min_elem=min(theta_l(k,:));
    max_elem=max(theta_l(k,:));
    columns_to_fill=find(theta_l(k,:)==0);
    filler=linspace(min_elem,max_elem,length(columns_to_fill)+2);
    theta_l(k,theta_l(k,:)==0)=filler(2:end-1);
end

% check that the nominal point theta_0 lies in the convex hull of each
% parameter combination of the given points
for k=1:length(kept_param_unique)
    params=kept_param_unique{k};
    if length(params)>=2 % for at least 2 dimensions
        checkhull=inhull(theta_0(params)',theta_l(params,:)');
    else % for 1 dimension
        if min(theta_l(params,:))<=theta_0(params) && theta_0(params)<=max(theta_l(params,:))
            checkhull=true;
        else
            checkhull=false;
        end
    end
    if checkhull==0
        warning(['theta_0 for parameters [',num2str(params),'] does not lie in the convex hull of the support points']);
    end
end

indata.theta_l=theta_l;

end

function best_combos=find_best_combo(kept_param_unique,unique_same_dim_positions,dim_unique)

permutations=cell(1,length(unique_same_dim_positions));
for k=1:length(unique_same_dim_positions) % loop for each parameter combination in a dimension e.g. [1 4], [2 4], [3 4]
    % find all possible permutations of each parameter combination
    % in given dimension e.g. {[1,4],[4,1]}, {[2,4],[4,2]},
    % {[3,4],[4,3]}
    permutations{k}=perms(kept_param_unique{unique_same_dim_positions(k)});
end

% find all possible combinations of permutations
% e.g.
% 1) [1,4;2,4;3,4]
% 2) [4,1;2,4;3,4]
% ...
rows=length(unique_same_dim_positions);
columns=dim_unique;
pages=size(permutations{1},1)^rows;
combinations=zeros(rows,columns,pages);
inputs=repmat({1:size(permutations{1},1)},1,rows);
outputs=allcomb(inputs{:});
scores=zeros(1,pages);
for k=1:pages
    for l=1:rows
        combinations(l,:,k)=permutations{l}(outputs(k,l),:);
    end
    
    % assign a score to each combination based on how many
    % same parameters are in the same column
    % e.g.
    % 1) [1,4;2,4;3,4] -> takes a score of 3 (all fours in the same column)
    % 2) [4,1;2,4;3,4] -> takes a score of 2 (only two fours in the same column)
    % ...
    for l=1:columns
        temp=accumarray(combinations(:,l,k),1);
        temp=temp(temp>1);
        scores(k)=scores(k)+sum(temp);
    end
end

% % choose the combination with the highest score
% [~,best_combo_index]=max(scores);
% best_combos=combinations(:,:,best_combo_index);

% choose the combinations with the (same) highest score
best_combos=combinations(:,:,scores==max(scores));

end

function theta_l_dim=make_theta_l_dim(theta_l_simplex,indata,kept_param_unique,unique_same_dim_positions,combo_cell)

n_theta=indata.n_theta;

theta_l_dim=cell(1,length(unique_same_dim_positions));

% apply smart rules to use as few support points as possible
for k=1:length(unique_same_dim_positions) % loop for each dimension e.g. 3D, 2D, 1D
    
    combo=combo_cell{k};
    
    % re-order rows of theta_l based on best_combo
    for l=1:length(unique_same_dim_positions{k})
        theta_l_simplex{unique_same_dim_positions{k}(l)}(combo(l,:),:)=theta_l_simplex{unique_same_dim_positions{k}(l)}(kept_param_unique{unique_same_dim_positions{k}(l)},:);
    end
    
    % "merge" together support points of each combination (as much
    % as possible)
    theta_l_dim{k}=repmat({zeros(n_theta,size(theta_l_simplex{unique_same_dim_positions{k}(1)},2))},1,length(unique_same_dim_positions{k})); % initialize
    %theta_l_dim{k}=repmat({zeros(n_theta,dim_unique(k)+1)},1,length(unique_same_dim_positions{k})); % initialize
    
    for l=1:length(unique_same_dim_positions{k})
        theta_l_index=unique_same_dim_positions{k}(l);
        
        theta_l_source=theta_l_simplex{theta_l_index};
        
        % find which parameters make the combination
        params=kept_param_unique{theta_l_index}; % best_combo(l,:);
        
        % check if existing support points contain any of the source points  
        source_points=theta_l_source(params,:);
        all_zero=false;
        m=1;
        n=1;
        %sum_same=0;
        while m<k && all_zero==false
            while n<=length(unique_same_dim_positions{m}) && all_zero==false
                destination_points=theta_l_dim{m}{n}(params,:);
                
%                 % find which rows contain non-zero elements
%                 nnz_flag=any(destination_points~=0,2);
            
                
                for op=1:size(source_points,2)
                    source_column=source_points(:,op);
                    qr=1;
                    point_found=false;
                    while qr<=size(destination_points,2) && point_found==false
                        destination_column=destination_points(:,qr);
                        nnz_flag=destination_column~=0;
                        
                        if source_column(nnz_flag)==destination_column(nnz_flag)
                            theta_l_dim{m}{n}(params,qr)=theta_l_source(params,op);
                            destination_points=theta_l_dim{m}{n}(params,:); % update destination points
                            
                            theta_l_source(params,op)=0; % zero means that these points already exist                            
                            point_found=true;
                        end
%                         [flag,column_index]=ismember(source_points(nnz_flag,op)',destination_points(nnz_flag,:)','rows');
%                         theta_l_dim{m}{n}(params,column_index(flag))=theta_l_source(params,op(flag));
%                         theta_l_source(params,op(flag))=0; % zero means that these points already exist
                        qr=qr+1;
                    end
                end
%                 %  % find the parameters that correspond to the non-zero rows
%                 %  params_to_check=params(nnz_flag);
%                 
%                 [same_points,point_inds]=ismember(source_points(nnz_flag,:)',destination_points(nnz_flag,:)','rows');
%                 
%                 theta_l_dim{m}{n}(params,point_inds(same_points))=theta_l_source(params,same_points);
%                 theta_l_source(params,same_points)=0; % zero means that these points already exist
%                 %sum_same=sum_same+sum(same_points);
                
                
                if all(theta_l_source==0,'all') %sum_same>=size(source_points,2)
                    all_zero=true;
                end
                n=n+1;
            end
            m=m+1;
        end
        
        % find the zero-columns of source and don't check these columns for
        % equality with destination columns
        points_to_check=any(theta_l_source~=0,1);
        
        % check if theta_l_source can be merged with existing
        % support points. If not, there is an empty space at
        % theta_l_temp{l}.
        m=1;
        merge_ok=false;
        while m<=l && merge_ok==false
            theta_l_dest=theta_l_dim{k}{m};
            
            % find the rows corresponding to parameters
            rows_to_check=theta_l_dest(params,:);
            
            % find which rows contain non-zero elements
            nnz_flag=any(rows_to_check~=0,2);
            
            % find the parameters that correspond to the non-zero rows
            params_to_check=params(nnz_flag);
            
            % if the non-zero rows are equal -> parameter combinations can be merged
            if isequal(theta_l_source(params_to_check,points_to_check),theta_l_dest(params_to_check,points_to_check))
                theta_l_dim{k}{m}(params,:)=theta_l_source(params,:);
                merge_ok=true;
            end
            m=m+1;
        end
    end
end



end

function theta_l=make_theta_l_merged(theta_l,n_theta,unique_same_dim_positions)

for k=2:length(unique_same_dim_positions) % loop for each dimension
    
    for l=1:length(unique_same_dim_positions{k}) % loop for each parameter combination in dimension
        theta_l_source=theta_l{k}{l};
        
        % find the non-zero rows of source = active parameters
        params=any(theta_l_source~=0,2);
        m=1;
        combination_checked=false;
        while m<=k && combination_checked==false && any(theta_l_source~=0,'all')
            n=1;
            while n<=length(unique_same_dim_positions{m}) && combination_checked==false
                theta_l_dest=theta_l{m}{n};
                
                % find the rows corresponding to parameters
                rows_to_check=theta_l_dest(params,:);
                
                % find how many zeros are in destination rows
                num_zeros=sum(rows_to_check==0,2);
                
                % move rows from source to destination if the number of
                % zeros in all rows of destination is at least equal to
                % the number of columns of source
                if all(num_zeros>=size(theta_l_source,2))
                    % find index of leftmost zero in destination rows
                    first_zero_index=nnz(rows_to_check(1,:))+1;
                    
                    theta_l{m}{n}(params,first_zero_index:first_zero_index+size(theta_l_source,2)-1)=theta_l_source(params,:);
                    theta_l{k}{l}(params,:)=0;
                    
                    combination_checked=true;
                end
                n=n+1;
            end
            m=m+1;
        end
    end
end
theta_l=[theta_l{:}]; % concatenate cells horizontally
theta_l=[theta_l{:}]; % concatenate matrices horizontally
theta_l(:,all(theta_l==zeros(n_theta,1),1))=[]; % remove zero columns

end

