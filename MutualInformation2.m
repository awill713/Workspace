function [ s2rInfoTemp MI ] = MutualInformation2( resp, err )

%MUTUALINFORMATION2 Summary of this function goes here
%   Calculates mutual information (MI) based on the responses and the
%   standard deviations associated with each response. MI is calculated
%   MI = sum(p(r|s) * log2(p(r|s) / p(r))). More info found in "Box 1" of:
%   https://www.nature.com/articles/nn1199_947.pdf

responses = resp;
stDev = err;

uniqueStimuli = length(responses);

%calculate overall probability of response based on Gaussian probability of
%response for each unique stimulus
lowest = min(responses - 2*stDev); %lowest response
highest = max(responses + 2*stDev); %highest response
increment = (highest-lowest)/1000;

%Variable that holds probability of each response
totalRespProb = zeros(2,length(lowest:increment:highest)); %row 1: response value, row 2: probability of that response value

totalRespProb(1,:) = lowest:increment:highest; %make row 1 the response values

%find probability of each response value (given the distriubtion of
%responses to all stimuli)
for r = 1:length(totalRespProb)
    
    tempR = 0; %temporary variable that sums the probability of response "r" for each unique stimulus
    
    for s = 1:uniqueStimuli
        %mean and std for Gaussian probability distribution
        avg = responses(s);
        sDev = stDev(s);
        
        %probability of response "r" given gaussian with center "avg" and std "sDev"
        dR = 1/(sDev*sqrt(2*pi)) * exp(-0.5*((totalRespProb(1,r)-avg)/sDev)^2);
        
        tempR = tempR + dR; %update tempR with dR
    end
    
    totalRespProb(2,r) = tempR / uniqueStimuli; %store tempR in the storage variable
end
%scale the overall response probability P(r) so it adds to 1
totalRespProb(2,:) = totalRespProb(2,:)/sum(totalRespProb(2,:));


%calculate the scalar that makes the response probability for each unique
%stimulus P(r|s) add to 1
probabilityScalar = zeros(1,uniqueStimuli); %scalar for each stimulus
% figure;
for s = 1:uniqueStimuli
    
    %Gaussian mean and std that the response probability is based on
    avg = responses(s);
    sDev = stDev(s);
    
    rx = lowest:increment:highest; %range of response probabilities r to sum over
    yx = zeros(1,length(rx)); %P(r|s) for each r
    
    %sweep through responses r and find P(r|s)
    for rr = 1:length(rx)
        yx(rr) = 1/(sDev*sqrt(2*pi)) * exp(-0.5*((rx(rr)-avg)/sDev)^2);
    end
    probScalar = 1/sum(yx); %scalar that makes sum of P(r|s)=1
    probabilityScalar(s) = probScalar;
%     yx = yx*probScalar;
%     plot(rx,yx);
%     hold on;
end

% figure;
% plot(totalRespProb(1,:),totalRespProb(2,:));


%calculate MI for each stimulus. MI = sum(p(r|s) * log2(p(r|s) / p(r)))
s2rInfoTemp = zeros(1,uniqueStimuli); %MI at each unique stimulus

for s = 1:uniqueStimuli
    tempS2Rinfo = 0; %variable that holds the incremental sum of MI of R|s
    
    %mean and std of Gaussian probability distrubution of responses
    theMean = responses(s);
    theSTD = stDev(s);
    scalar = probabilityScalar(s); %scalar that makes the gaussian probability distribution sum to 1
    
    %MI = p(r|s) * log2(p(r|s) / p(r))
    for r = 1:length(totalRespProb)
        
        theRgivenS = scalar * 1/(theSTD*sqrt(2*pi)) * exp(-0.5*((totalRespProb(1,r)-theMean)/theSTD)^2);
        
        dMI = theRgivenS * log2(theRgivenS / totalRespProb(2,r));
        tempS2Rinfo = tempS2Rinfo + dMI;
    end
    s2rInfoTemp(s) = tempS2Rinfo;
end

MI = sum(s2rInfoTemp); %overall MI is sum of MI for each unique stimulus

% figure;
% plot(s2rInfoTemp);

end

