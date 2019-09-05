%%
% Perform network analysis on spontaneous and tone-evoked data
spontNetworkAnalysis
toneEvokedNetworkAnalysis

% Are driven units statistically more likely to be in any given cluster
% than random, using the clustering from spontaneous data (netanal)
drivenUnitsFunctionalModules

% Are tuned units statistically more likely to be in any given cluster than
% random, using the clustering from tone-evoked data (netanal_tone)
tunedUnitsFunctionalModules
% using spont clusters (netanal)
tunedUnitsFunctionalModules_spont

% Are units with the same best frequencies more likely than chance to be in
% the same clusters, using spontaneous clusters (netanal)
bestFreqFunctionalModules
% using tone-evoked clusters (netanal_tone)
bestFreqFunctionalModules_toneEvoked

% Do units in the same cluster have a higher correlation in their frequency
% responses than other cluster/random?
tuningCurveCorrelation_modules % NOT FINISHED

% save ordered correlations and image of clustered units on imaging plane
saveNetworkAnalFigs

% Plot tuning curves of 'good cells' in each cluster

%% 28 Nov 2017 update - Rick has sent new cluster labels. I will look into them to see if there are any interesting things...

drivenUnitModules





