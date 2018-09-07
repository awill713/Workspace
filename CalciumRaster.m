
%you need the calcium file from "...processedData"

ystep = 10000;

figure;
hold on;
for i = 1:length(tunedNeurons)
    offset = ystep*i;
    plot(calcium.npilSubTraces(tunedNeurons(i),:)+offset);
end
totalTime = length(calcium.npilSubTraces)/30/60; %30 fps, 60 sec per min
xticks((3:3:totalTime+3)*30*60);
xticklabels({'3','6','9','12','15','18'});
xlabel('Time (min)');
yticks(5*ystep:5*ystep:offset+5*ystep);
yticklabels({'5','10','15','20','25'});
ylabel('Neuron number');