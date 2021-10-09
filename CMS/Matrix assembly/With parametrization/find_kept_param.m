function kept_param=find_kept_param(dofdata,indata,interface)
    adj_I=dofdata.adj_I;
    bound_S=dofdata.bound_S;
    group_l=indata.group_l;
    N_S=indata.N_S;
    S_0=indata.S_0;
    S_j=indata.S_j;
    boundaries=vertcat(adj_I{group_l{interface}});
    kept_param=[];
    for k=1:N_S 
        if ~isempty(intersect(boundaries,bound_S{k}))
            param=[];
            if ismember(k,S_0)
                %param=0;
                % if interface is adjacent to a group which does not depend to any model parameter,
                % do not add any parameter (zero has no meaning and breaks the code)
                param=[]; 
            else
                l=0;
                while isempty(param)
                    l=l+1;
                    if ismember(k,S_j{l})
                        param=l;
                    end
                end
            end
            kept_param=[kept_param param];
        end
    end
    kept_param=unique(kept_param);
end