function indata=error_check(dofdata,indata)

n_id_S=indata.n_id_S;
n_IR_l=indata.n_IR_l;
N_S=indata.N_S;
reduction_I=indata.reduction_I;
method_eigf_S=indata.eigf.group.method;
method_eigf_I=indata.eigf.interface.method;

N_I=dofdata.N_I;

if length(n_id_S)<N_S
    error('Check dimension of n_id_S.');
elseif length(n_id_S)>N_S
    warning('Extra values of n_id_S will be ignored.');
    indata.n_id_S(N_S+1:end)=[]; % extra values must be deleted or else there will be error in matrix assembly
    indata.n_id=sum(indata.n_id_S);
end

if method_eigf_S==1 % perform additional check only if cutoff freuqencies are used
    target=indata.eigf.group.target;
    max=indata.eigf.group.max;
    step=indata.eigf.group.step;
    init=indata.eigf.group.init;
    if length(target)<N_S || length(max)<N_S || length(step)<N_S || length(init)<N_S % one of these contains less values than # of component groups
        error('Check dimension of target_eigf_S, max_eigf_S, step_eigf_S and init_eigf_S.');
    end
end

if reduction_I==3 || reduction_I==6 % local interface reduction is used
    if length(n_IR_l)<N_I
        error('Check dimension of n_IR_l.');
    elseif length(n_IR_l)>N_I
        warning('Extra values of n_IR_l will be ignored.');
        indata.n_IR_l(N_I+1:end)=[]; % extra values must be deleted or else there will be error in matrix assembly
        indata.n_IRL=sum(indata.n_IR_l);
    end    
end

if reduction_I~=1 && reduction_I~=4 % interface reduction is used (global or local)
    if method_eigf_I==1 % perform additional check only if cutoff freuqencies are used
        target=indata.eigf.interface.target;
        max=indata.eigf.interface.max;
        step=indata.eigf.interface.step;
        init=indata.eigf.interface.init;
        if length(target)<N_I || length(max)<N_I || length(step)<N_I || length(init)<N_I % one of these contains less values than # of interfaces
            error('Check dimension of target_eigf_I, max_eigf_I, step_eigf_I and init_eigf_I.');
        end
    end
end

indata.n_D=indata.n_id+dofdata.n_I;
indata.n_DI=indata.n_id+indata.n_IR;
indata.n_DIL=indata.n_id+indata.n_IRL;

end