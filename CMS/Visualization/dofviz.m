function dofviz(moddata,dofdata,indata,viz)
xmsh=moddata.xmsh;
dim=moddata.dim;
N_S=indata.N_S;
i_dofs_S=dofdata.i_dofs_S;
b_dofs_S=dofdata.b_dofs_S;
I_dofs_l=dofdata.I_dofs_l;
N_I=dofdata.N_I;

if string(viz)=="all" % display dofs for each group only if viz is set to "all"
    for k=1:N_S
        figure;
        hold on
        if dim==2
            coords=findcoords(xmsh,i_dofs_S{k});
            scatter(coords(1,:),coords(2,:),'b.');
            coords=findcoords(xmsh,b_dofs_S{k});
            scatter(coords(1,:),coords(2,:),'r.');        
        else
            coords=findcoords(xmsh,i_dofs_S{k});
            scatter3(coords(1,:),coords(2,:),coords(3,:),'b.');
            coords=findcoords(xmsh,b_dofs_S{k});
            scatter3(coords(1,:),coords(2,:),coords(3,:),'r.');
            view(3)
        end
        axis equal
        grid on
        title(['Group ',num2str(k)]);
        legend('Internal dofs', 'Boundary dofs');
    end
end

% always display interface dofs
figure;
hold on
leg={}; % define in case there are no boundaries (a single domain) 
for k=1:N_I
    coords=findcoords(xmsh,I_dofs_l{k});  
    if dim==2
        scatter(coords(1,:),coords(2,:),'.');
    else
        scatter3(coords(1,:),coords(2,:),coords(3,:),'.');
        view(3)
    end    
    leg{k}=['Interface ',num2str(k)];
end
axis equal
grid on
title('Interfaces');
legend(leg);

end

function coords=findcoords(xmsh,dofs_S)
    %[linindex,~]=find(xmsh.nodes.dofs(:)==dofs_S');
    [~,linindex]=intersect(xmsh.nodes.dofs(:),dofs_S);
    [~,col]=ind2sub(size(xmsh.nodes.dofs),linindex);
    col=unique(col);
    coords=xmsh.nodes.coords(:,col);
end