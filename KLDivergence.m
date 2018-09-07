function [ kld ] = KLDivergence( x, y )
%KLDIVERGENCE Summary of this function goes here
%   Detailed explanation goes here

bins = 1000;

low = min([min(x) min(y)]);
high = max([max(x) max(y)]);
range = linspace(low,high,bins+1);

[xMean, xSigma] = normfit(x);
xNorm = normpdf(range,xMean,xSigma);
xProb = xNorm./sum(xNorm);
[yMean, ySigma] = normfit(y);
yNorm = normpdf(range,yMean,ySigma);
yProb = yNorm./sum(yNorm);
% figure;hold on;
% plot(range,xNorm);
% plot(range,yNorm);

kld = sum(-xProb.*log2(yProb./xProb));

% xHist = histc(x,linspace(low,high,bins+1));
% xProb = xHist./sum(xHist);
% yHist = histc(y,linspace(low,high,bins+1));
% yProb = yHist./sum(yHist);
% 
% figure; hold on;
% plot(xProb);
% plot(yProb);
% 
% kld = 0;
% for b = 1:bins
%     if xProb(b) == 0
%         delta_kld = 0;
%     else
%         delta_kld = -1 * xProb(b) * log2(yProb(b) / xProb(b));
%     end
% %     delta_kld = -1 * xProb(b) * log2(yProb(b) / (xProb(b) + eps));
%     kld = kld + delta_kld;
% end
% 
% kld2 = sum(-xProb.*log2(yProb./xProb+eps));

end

