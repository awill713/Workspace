function [ s2rInfoTemp MI ] = MutualInformation2( resp, err )
%MUTUALINFORMATION2 Summary of this function goes here
%   Detailed explanation goes here

responses = resp;
stDev = err;

uniqueStimuli = length(responses);
mutInf = zeros(length(responses));


%calculate overall probability of response
lowest = min(responses - 2*stDev);
highest = max(responses + 2*stDev);
increment = (highest-lowest)/1000;
totalRespProb = zeros(2,length(lowest:increment:highest)); %response value and probability of that response value

totalRespProb(1,:) = lowest:increment:highest;
for r = 1:length(totalRespProb)
    tempR = 0;
    
    for s = 1:uniqueStimuli
        avg = responses(s);
        sDev = stDev(s);
        dR = 1/(sDev*sqrt(2*pi)) * exp(-0.5*((totalRespProb(1,r)-avg)/sDev)^2);
        
        tempR = tempR + dR;
    end
    totalRespProb(2,r) = tempR / uniqueStimuli;
end

%scale the probabilities for each stimulus such that they add to 1
% figure;
for s = 1:uniqueStimuli
    avg = responses(s);
    sDev = stDev(s);
    rx = lowest:increment:highest;
    yx = zeros(1,length(rx));
    for rr = 1:length(rx)
        yx(rr) = 1/(sDev*sqrt(2*pi)) * exp(-0.5*((rx(rr)-avg)/sDev)^2);
    end
    probScalar = 1/sum(yx);
    probabilityScalar(s) = probScalar;
    yx = yx*probScalar;
%     plot(rx,yx);
%     hold on;
    sum(yx);
end
% legend('1','2','3','4','5','6','7','8','9','10');

totalRespProb(2,:) = totalRespProb(2,:)/sum(totalRespProb(2,:));
% figure;
% plot(totalRespProb(1,:),totalRespProb(2,:));

% totalRespProb(1,:)
% totalRespProb(2,:)

%calculate MI for each stimulus
s2rInfoTemp = zeros(1,uniqueStimuli);

for s = 1:uniqueStimuli
    tempS2Rinfo = 0;
    theMean = responses(s);
        theSTD = stDev(s);
        scalar = probabilityScalar(s);
    for r = 1:length(totalRespProb)
        
        theRgivenS = scalar * 1/(theSTD*sqrt(2*pi)) * exp(-0.5*((totalRespProb(1,r)-theMean)/theSTD)^2);
        
        dMI = theRgivenS * log2(theRgivenS / totalRespProb(2,r));
        tempS2Rinfo = tempS2Rinfo + dMI;
    end
    s2rInfoTemp(s) = tempS2Rinfo;
end

MI = sum(s2rInfoTemp);

% figure;
% plot(s2rInfoTemp);

end

