

repeats = 1000;
intensity = 0.1:0.1:10;
output = zeros(1,length(intensity));
dev = zeros(1,length(intensity));
var = zeros(1,length(intensity));
fano = zeros(1,length(intensity));
cv = zeros(1,length(intensity));

for i = 1:length(intensity)
    trials = poissrnd(intensity(i),[1 repeats]);
    stdd = std(trials);
    avg = mean(trials);
    output(i) = avg;
    dev(i) = stdd;
    var(i) = stdd^2;
    fano(i) = stdd^2/avg;
    cv(i) = stdd/avg;
end

figure;hold on;
hA = area(intensity,[output-dev; 2*dev]');
hA(1).FaceAlpha = 0;
hA(1).EdgeColor = [1 1 1];
hA(2).FaceColor = [0 0 0];
hA(2).FaceAlpha = 0.5;
hA(2).EdgeColor = [1 1 1];
plot(intensity,output,'Color',[0 0 0]);
xlabel('Poisson \lambda');
ylabel('Firing rate (Hz)');
