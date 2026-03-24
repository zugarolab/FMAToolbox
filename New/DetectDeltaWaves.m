function deltaWaves = DetectDeltaWaves(channel,spikes,save)

%DetectDeltaWaves - Detect cortical delta waves from LFP and refine selection using local spiking activity.
%
%    Detects delta waves (1–6 Hz) from LFP signals by:
%      1) Extracting LFP from the specified channel (GetLFP)
%      2) Cleaning the signal to remove artifacts (CleanLFP)
%      3) Detecting candidate delta waves (FindDeltaWaves)
%      4) Applying an amplitude-based selection using peak–trough difference
%
%    If spike data are provided, the function enables manual refinement
%    of the threshold by visualizing peri-event time histograms (PETH)
%    of spike activity aligned to delta wave peaks.
%
%  USAGE
%
%    deltaWaves = DetectDeltaWaves(channel,spikes,save)
%
%    channel        LFP channel index used to extract the signal.
%    spikes         optional spike times. If provided, allows interactive 
%                   threshold refinement based on spike activity around 
%                   delta waves.
%    save           optional boolean (default = true). If true, results are
%                   saved to disk as a .deltaWaves.events.mat file.
%
%  OUTPUT
%
%    deltaWaves     structure containing detected delta wave properties:
%                   - timestamps          [start end] times of each wave
%                   - peaks               peak times of each wave
%                   - peakNormedPower     z-scored peak amplitudes
%                   - troughValue         z-scored trough amplitudes
%                   - detectorName        description of detection method
%                   - badIntervals        intervals removed during cleaning
%
%  NOTE
%
%    The peak–trough threshold is session-dependent and should ideally be
%    adjusted manually when spike data are available. Default threshold is 4.5
%    (z-scored units).
%
%  SEE
%
%    See also GetLFP, CleanLFP, FindDeltaWaves, PETH.

basepath = pwd;
[~,basename] = fileparts(basepath);
if ~exist('save','var'), save = true; end

lfp = GetLFP(channel);
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 Inf],'manual',true);
deltas0 = FindDeltaWaves(clean);

threshold = 4.5;
if exist('spikes','var') && ~isempty(spikes)
    nQuantiles = 100; durations = [-1 1]; nBins = 201; smooth = [0 1];
    figure;
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

