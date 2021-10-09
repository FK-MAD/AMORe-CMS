function correspondence_I=geomviz(model,dofdata,indata,groupviz,renumber_I)

figure;

mphgeom(model,'geom1','domainlabels','on','domainlabelscolor','r',...
    'facelabels','on','facelabelscolor','g','facealpha',0);
 
% find text objects with green color (corresponding to face labels)
g_array=findobj('Type','text','Color','g'); 

g_txt=str2double({g_array(:).String});

% find the text objects to hide (those that are not interfaces)
hiden=setdiff(g_txt,dofdata.bound_I,'stable');
[~,g_txt_index]=intersect(g_txt,hiden,'stable');

set(g_array(g_txt_index),'visible',0);

red_text='\fontsize{8} {\color{red}red text} \rightarrow domains';
green_text='{\color{green}green text} \rightarrow interfaces (original numbering)';



if string(groupviz)=="groups"
    
    % find text objects with red color (corresponding to face labels)
    red_array=findobj('Type','text','Color','r'); 
    set(red_array,'visible',0); % hide all domain numbers
    r_txt=str2double({red_array(:).String});
    
    N_S=indata.N_S;
    group_S=indata.group_S;
    for k=1:N_S
        [~,r_txt_index]=intersect(r_txt,group_S{k},'stable');
        pos=vertcat(red_array(r_txt_index).Position);
        pos=sum(pos,1)/size(pos,1);
        text(pos(1),pos(2),pos(3),num2str(k),'color','r'); % number of group in the centroid of all domains which belong to the group
    end
    red_text='\fontsize{8} {\color{red}red text} \rightarrow groups';
end

% define in case there is no renumbering
correspondence_I=[];
if string(renumber_I)=="renumber_I"
    g_array=findobj('Type','text','Color','g','Visible','on');
    [correspondence_I,sort_index]=sort(str2double({g_array(:).String}));
    set(g_array(sort_index),{'String'},num2cell(1:length(g_array))');
    green_text='{\color{green}green text} \rightarrow interfaces (renumbered)';
end

title({'Geometry (geom1)';red_text;green_text});


end