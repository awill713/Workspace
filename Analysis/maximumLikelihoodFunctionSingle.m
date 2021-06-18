function [outputDir] = maximumLikelihoodFunctionSingle(unitStats,testProbe)
%MAXIMUMLIKELIHOODFUNCTIONSINGLE Summary of this function goes here
%   Detailed explanation goes here

likelihoods = zeros(1,size(unitStats,1));
for like = 1:length(likelihoods)
    post = pdf(unitStats(like),testProbe);
    
    likelihoods(like) = log(post);
end
% likelihoods
[outputProb, outputDir] = max(likelihoods);

end

