clear clusterChannelAndDepth 

mouseID = 'AW118';
session = 'Session3';
folder = '2020-02-21_19-54-33';

%where to save the unit depth information
outDir = fullfile('D:\KiloSort',mouseID,session,folder);

%load the needed files
spike_clusters = double(readNPY([outDir filesep 'recording01\spike_clusters.npy'])); %which unit (after manually sorted) each spike belongs to
spike_templates = double(readNPY([outDir filesep 'recording01\spike_templates.npy'])); %which template (Kilosort-assigned units before manually sorting) each spike belongs to
cluster_group = importdata([outDir filesep 'recording01\cluster_group.tsv'],'\t'); %final unit numbers and assignments (noise, good, multi), first row is just info
templates = double(readNPY([outDir filesep 'recording01\templates.npy'])); %the average spike waveform of each template (Kilsort-assigned units before manually sorting) across all channels
channel_positions = double(readNPY([outDir filesep 'recording01\channel_positions.npy']));  %depth info of each channel


%get the KiloSort/Phy unit numbers of the final manually-sorted units
labels = vertcat(cluster_group(2:end));

%for each unit
for c = 1:length(labels)
    
    %get the final KiloSort/Phy unit number of this unit, after all manual
    %sorting/splitting/merging is done
    str = strsplit(labels{c}); 
    ID = str2num(str{1});
    
    %get the templates (Kilosort-assigned, pre manually sorted units) that
    %this unit came from. Usually just a 1:1. But one unit could come from
    %multiple templates if you manually merged the units when sorting
    templatesOfInterest = unique(spike_templates(spike_clusters==ID)) + 1; % +1 because matlab is 1-indexed, whereas Phy is 0-indexed
    
    %the template amplitude on each channel is the peak-to-trough (max -
    %min) on each channel. Average across the templates that comprise this
    %unit in order to get the unit's amplitude across all channels 
    amp = mean(max(templates(templatesOfInterest,:,:),[],2) - min(templates(templatesOfInterest,:,:),[],2),1);
    
    %the assigned channel is where the amplitude is the largest
    clusterChannelAndDepth(c) = find(amp==max(amp));
end

%the depth is assigned according to the probe map
clusterChannelAndDepth(2,:) = channel_positions(clusterChannelAndDepth(1,:),2);

save([outDir filesep 'unitDepths'],'clusterChannelAndDepth');

