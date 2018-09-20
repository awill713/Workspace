function [ rampedSignal ] = applyRamp( train, dur)
%APPLYRAMP Summary of this function goes here
%   Detailed explanation goes here

rampedSignal = train;
rampDur = dur; %duration in samples (NOT milliseconds)

for r = 1:rampDur
    rampCoef = pi / (2*rampDur);
    rampedSignal(r) = sin(rampCoef * r) * rampedSignal(r);
    
    tailIndex = length(rampedSignal) - rampDur + r;
    rampedSignal(tailIndex) = sin(rampCoef*(r-1) + pi/2) * rampedSignal(tailIndex);
end

end

