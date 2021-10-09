function [PHI,LAMBDA,kept]=keptmodes(K,M,method,kept,target,max,step,init)

if method==0
    fprintf(['        Calculating lowest ',num2str(kept),' modes...\n']);
    [PHI,LAMBDA]=eigs(K,M,kept,'smallestabs');
else
    n_calc=init;
    n_smaller=0;
    LAMBDA=[];
    k=0;
    while n_smaller==length(diag(LAMBDA))
        if n_calc>max
            error(['Reached maximum allowed number of kept modes (',num2str(max),').',...
                ' Consider decreasing the step and/or increasing the maximum number of kept modes.']);
        end
        k=k+1;
        fprintf(['        Attempt #',num2str(k),' at finding lowest modes up to cutoff frequency ',num2str(target),' Hz...\n']);
        fprintf(['        Calculating lowest ',num2str(n_calc),' modes...\n']);
        [PHI,LAMBDA]=eigs(K,M,n_calc,'smallestabs');
        n_calc=n_calc+step;
        n_smaller=sum(diag(LAMBDA)<(target*2*pi)^2);
    end
    kept=n_smaller+1;
    fprintf(['        Cutoff frequency reached. Kept the lowest ',num2str(kept),' modes...\n']);
    PHI=PHI(:,1:kept);
    LAMBDA=LAMBDA(1:kept,1:kept);
end

end