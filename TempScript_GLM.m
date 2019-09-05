
clear

%Independent variable:    1     2   3
%Dependent variable mean: 100   0   50
repeats = 100000;
elements = repeats * 4;

% independent = datasample([1 2],elements);
% independent = randi(2,[2 elements])-1;

order = repmat([1 2 3 4],[1 repeats]);
order = order(randperm(length(order)));

independent = zeros(2,elements);
for e = 1:elements
    if order(e)==1
        independent(:,e) = [1;0];
    elseif order(e)==2
        independent(:,e) = [0;1];
    elseif order(e)==3
        independent(:,e) = [1;1];
    elseif order(e)==4
        independent(:,e) = [0;0];
    end
end

for e = 1:elements
%     dependent(e) = 100*independent(1,e) + 10*independent(2,e) + 50*independent(3,e) + randn(1)*20;
    if independent(1,e) && independent(2,e)
        dependent(e) = 0;
    elseif independent(1,e)
        dependent(e) = 1;
    elseif independent(2,e)
        dependent(e) = 1;
    else
        dependent(e) = 0;
    end
    dependent(e) = dependent(e) + randn(1)*0.001;
end

[beta sigma resid] = mvregress(independent',dependent');
beta

glm = fitglm(independent',dependent','interactions','distr','normal');

[p table stats] = anova2([dependent(find(order==4))' dependent(find(order==1))';dependent(find(order==2))' dependent(find(order==3))'],repeats,'off');