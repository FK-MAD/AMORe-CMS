function MACviz(exact,approximate,num)

figure
for k=1:num
    for l=1:k
        ex=exact(:,k);
        apr=approximate(:,l);
        MAC(l,k)=(ex'*apr*apr'*ex)/(ex'*ex*apr'*apr);
    end
end

MAC=MAC+triu(MAC,1)';
MAC=rot90(MAC);
h=bar3(MAC);

view(45,45);
c=colorbar;
c.Label.String='MAC Value';
colormap jet

xlim([0 num+1]);
xticks(1:num);
xlabel('Mode Number');

ylim([0 num+1]);
yticks(1:num);
ylabel('Mode Number');
yticks_reverse = arrayfun(@num2str,num:-1:1,'uni',false);
yticklabels(yticks_reverse)

% code below is from:
% https://www.mathworks.com/matlabcentral/answers/98236-how-can-i-color-bars-to-correspond-to-their-heights-when-using-bar3#answer_107586
numBars = size(MAC,1);
numSets = size(MAC,2);
for i = 1:numSets
    zdata = ones(6*numBars,4);
    k = 1;
    for j = 0:6:(6*numBars-6)
        zdata(j+1:j+6,:) = MAC(k,i);
        k = k+1;
    end
    set(h(i),'Cdata',zdata)
end

end