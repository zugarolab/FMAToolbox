function deltaWaves = DetectDeltaWaves(channel,spikes,save)

basepath = pwd;
[~,basename] = fileparts(basepath);
if ~exist('save','var'), save = true; end

lfp = GetLFP(channel);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 Inf],'manual',true);
deltas0 = FindDeltaWaves(clean);

threshold = 4.5;
if exist('spikes','var') && ~isempty(spikes)
    figure;
    nQuantiles = 100; durations = [-1 1]; nBins = 201;
    strength = deltas0(:,5)-deltas0(:,6); 
    [h,ht] = PETH(spikes(:,1),deltas0(:,2),'durations',durations,'nBins',nBins);
    
    clf;
    p = mean(strength<threshold);
    PlotColorMap(Smooth(Shrink(sortby(h,strength),floor(size(strength,1)/nQuantiles),1),smooth),'x',ht);
    PlotHVLines(ceil(p*nQuantiles)+0.5,'h','w--');
    disp('Change the threshold manually');
    
    keyboard;
end

deltas = deltas0(deltas0(:,5)-deltas0(:,6)>threshold,:); % these thresholds should be manually refined for each session
deltaWaves.timestamps = deltas(:,[1 3]); deltaWaves.peaks = deltas(:,2);  
deltaWaves.peakNormedPower = deltas(:,5); 
deltaWaves.detectorName = ['channel ' num2str(channel) '(+1), CleanLFP, FindDeltaWaves, peak-trough> ' num2str(threshold)];
deltaWaves.troughValue = deltas(:,6); deltaWaves.badIntervals = badIntervals;

save(fullfile(basepath,[basename '.deltaWaves.events.mat']),'deltaWaves');
