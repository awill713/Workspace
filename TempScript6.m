
values = linspace(1,9,100);
curve = normpdf(values,5,2);

for i = 1:length(curve)
    temp = zeros(1,10);
    for j = 1:length(temp)
        temp(j) = curve(i);
    end
    n4r(i) = mean(temp);
    n4e(i) = std(temp);
    
    temp = zeros(1,10);
    for j = 1:length(temp)
        temp(j) = curve(i)*2;
    end
    n5r(i) = mean(temp);
    n5e(i) = std(temp);
    
    temp = zeros(1,10);
    for j = 1:length(temp)
        temp(j) = curve(i)*1.05;
    end
    n6r(i) = mean(temp);
    n6e(i) = std(temp);
end
figure;hold on
plot(n4r,'LineWidth',4,'Color','k');
plot(n5r,'LineWidth',4,'Color','b');
plot(n6r,'LineWidth',4,'Color','r');