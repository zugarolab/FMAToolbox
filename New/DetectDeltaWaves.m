function deltaWaves = DetectDeltaWaves(channel,spikes,doSave,opt)

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
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'restrict'    intervals (s) to restrict analysis to
%     'session'     path to session .xml file, delta waves file will be saved
%                   in its session folder (default: current location, assumed
%                   to be a session folder)
%     'fast'        if 'on', 'spikes' is assumed to be time sorted, increasing
%                   speed (defaul: 'off', ignored if 'spikes' is empty)
%     'check'       number of randomly chosen examples to display (default: 0)
%    =========================================================================
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

% Copyright (C) 2017-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    channel
    spikes = []
    doSave {mustBeGeneralLogical} = true
    opt.restrict (:,2) {mustBeNumeric} = []
    opt.session (1,1) string = ""
    opt.fast {mustBeGeneralLogical} = 'off'
    opt.check (1,1) {mustBeInteger,mustBeNonnegative} = 0
end

doSave = GeneralLogical(doSave);

% default values
if opt.session == ""
    basepath = pwd;
    [~,basename] = fileparts(basepath);
else
    [basepath,basename] = fileparts(opt.session);
end

lfp = GetLFP(channel,'restrict',opt.restrict);
% remove artifacts
[clean,~,badIntervals] = CleanLFP(lfp,'thresholds',[6 Inf],'manual',true);
% putative delta waves
deltas0 = FindDeltaWaves(clean);

strength = deltas0(:,5) - deltas0(:,6);
threshold = 4.5;
if ~isempty(spikes)
    % peth parameters
    durations = [-1 1];
    nBins = 201;
    smooth = [0 1];
    nQuantiles = 100;

    figure;
    if any(diff(deltas0(:,2))<= 0), disp("Unexpected time ordering of putative deltas waves, option 'fast' should not be used"), end
    [h,ht] = PETH(spikes(:,1),deltas0(:,2),'durations',durations,'nBins',nBins,'fast',opt.fast); % MUA vs Delta peak time peth
    clf;

    PlotColorMap(Smooth(Shrink(sortby(h,strength),floor(size(strength,1)/nQuantiles),1),smooth),'x',ht);
    plotThreshold = @(x) yline(ceil(mean(strength<x)*nQuantiles)+0.5,'k--'); % utility to see new thresholds
    plotThreshold(threshold)
    % give control to user
    fprintf(1,'Change "threshold" manually, then type:\n  "plotThreshold(threshold)" to visualize it,\n  "dbcont" to continue.\n');
    keyboard

end

% this threshold should be manually refined for each session
deltas = deltas0(strength>threshold,:);
deltaWaves.timestamps = deltas(:,[1 3]);
deltaWaves.peaks = deltas(:,2);  
deltaWaves.peakNormedPower = deltas(:,5); 
deltaWaves.detectorName = ['channel ' num2str(channel) ' (+1), CleanLFP, FindDeltaWaves, peak-trough > ' num2str(threshold)];
deltaWaves.troughValue = deltas(:,6);
deltaWaves.badIntervals = badIntervals;

% save
if doSave
    save(fullfile(basepath,basename+'.deltaWaves.events.mat'),'deltaWaves');
end

% see examples
if opt.check
    % choose random examples
    idx = randperm(numel(deltaWaves.peaks),opt.check);
    intervals = deltaWaves.peaks(idx) + [-.5,.5];
    figure
    for i = 1 : opt.check
        clf, hold on
        title(sprintf('%s Example %d/%d - l-click: next, r-click: exit',replace(basename,'_','\_'),i,opt.check));
        xlabel('time (s)'),
        ylabel('LFP')
        PlotXY(Restrict(lfp,intervals(i,:)))
        xline(deltaWaves.timestamps(idx(i),1),'Color','k','DisplayName','Start');
        xline(deltaWaves.peaks(idx(i)),'Color','g','DisplayName','Peak');
        xline(deltaWaves.timestamps(idx(i),2),'Color','r','DisplayName','Stop');
        xlim(intervals(i,:))
        [~,~,button] = ginput(1);
        if button == 3 % right click
            break
        end
    end

end