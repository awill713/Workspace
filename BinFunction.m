function [ binnedTrain ] = BinFunction(func, spikeTrain)
%BINFUNCTION Summary of this function goes here
%   Detailed explanation goes here

if func==2 %sliding bins, not continuous
    
    binnedTrain = zeros(1,length(spikeTrain));
    
    binLength=100;
    
    for i=binLength:length(spikeTrain)
        
        binnedTrain(i) = sum(spikeTrain(i-binLength+1:i))/binLength;
        
    end
    
else
    
    alpha = 100; %used to be 100
    
    for i=2:length(spikeTrain)
        
        intSum = 0;
        
        ti = i*0.001;
        
        for j=1:i-1
            
            tj = j*0.001;
            
            temp = 0.001*(alpha*alpha*tj*exp(-alpha*tj))*spikeTrain(i-j);
            
            if temp<0
                temp=0;
            end
            intSum = intSum + temp;
            
        end
        
        binnedTrain(i) = intSum;
        
    end
end

