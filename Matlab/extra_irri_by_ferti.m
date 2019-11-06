function extra_irri_by_ferti(fsol2, strategy)

theta_table = dlmread('theta_table.txt')
if strcmp (strategy, 'Fert')
for i=1:size(theta_table,1)
    index = find(fsol2(:,1)==theta_table(i,1))
    if isempty(index)==0
    theta_table(i,4) = theta_table(i,4)+ 0.2;
    end
end
end

dlmwrite('theta_table.txt', theta_table)
     
     