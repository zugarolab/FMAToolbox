function [clean,bad,badIntervals] = CleanLFP(lfp,opt)
% CleanLFP - Remove artefacts from LFP trace
%
% Returns the lfp cleaned from two kinds of artefacts: big artefacts in the
% lfp signal (surpassing a threshold in z-units) and very fast fluctuations,
% in which the signal derivative suprasses a threshold in z-units).
%
%  USAGE
%
%    [clean,bad,badIntervals] = CleanLFP(lfp,<options>)
%
%    data               lfp to clean in [timestamps values] format, where
%                       timestamps are expected to be in secods.
%    <options>          optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties        Values
%    -------------------------------------------------------------------------
%     'thresholds'      the thresholds for the two signal and derivative
%                       artefacts in z-units (default = [5 10])
%     'aroundArtefact'  time around artefact boundaries that should be removed,
%                       in seconds (default = [2 0.1])
%     'manual'          choose your thresholds manually in debug mode
%    =========================================================================
%
%  OUTPUT
%
%    clean          the lfp in which the rows inside artefacts have been deleted
%    bad            a logical vector indicating which rows of the original lfp file
%                   were flagged as within the artefact and should be removed
%    badIntervals   the [start stop] intervals in which the signal or derivative
%                   surpassed its respective threshold

% Copyright (C) 2017-2022 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
  lfp
  opt.thresholds (1,2) {mustBeNumeric} = [5,10]
  opt.aroundArtefact (1,2) {mustBeNumeric} = [0.5,0.1]
  opt.manual {mustBeGeneralLogical} = false
end

opt.manual = GeneralLogical(opt.manual);

threshold1 = opt.thresholds(1); % in sigmas deviating from the mean
aroundArtefact1 = opt.aroundArtefact(1); % 2, Big and long artefacts

threshold2 = opt.thresholds(2); % for derivative of z-scored signal
aroundArtefact2 = opt.aroundArtefact(2); % 0.1 Very fast fluctuations (short time scale)

t = lfp(:,1);
values = lfp(:,2);
z = zscore(values);
d = [diff(z);0];
bad = false(size(values));

if opt.manual
    figure, hold on
    plot(t,z);
    xlabel('time (s)'); ylabel('lfp signal (z-units)');
    plotThreshold = @(x) yline(x*[-1,1],'r--'); % utility to see new thresholds
    plotThreshold(threshold1);
    legend('signal','threshold1');
    fprintf(1,'Change "threshold1" manually, then type:\n  "plotThreshold(threshold1)" to visualize it,\n  "dbcont" to continue.\n');
    keyboard

    % Optionally, take a look at the threshold for the derivative
    if false
        clf, hold on;
        plot(t,d); plotThreshold(threshold2); ylabel('lfp derivative (z-units)'); legend('derivative','threshold2');
        fprintf(1,'Change "threshold2" manually, then type:\n  "plotThreshold(threshold2)" to visualize it,\n  "dbcont" to continue.\n');
    end
end

% Detect large global artefacts (1)
artefactInterval = t(FindInterval(abs(z) > threshold1));
if numel(artefactInterval) == 2
    artefactInterval = artefactInterval(:)';
end
if ~isempty(artefactInterval)
    artefactInterval = ConsolidateIntervals([artefactInterval(:, 1) - aroundArtefact1, artefactInterval(:, 2) + aroundArtefact1]);
    bad = InIntervals(t, artefactInterval);
else
    artefactInterval = zeros(0, 2);
end

% Find noise using the derivative of the zscored signal (2)
noisyInterval = t(FindInterval(abs(d) > threshold2));
if numel(noisyInterval) == 2
    noisyInterval = noisyInterval(:)';
end
if ~isempty(noisyInterval)
    noisyInterval = ConsolidateIntervals([noisyInterval(:, 1) - aroundArtefact2, noisyInterval(:, 2) + aroundArtefact2]);
    bad = bad | InIntervals(t, noisyInterval);
else
    noisyInterval = zeros(0, 2);
end

% Substitute noisy signal with interpolated signal as if artefact did not exist
values(bad) = interp1(t(~bad), values(~bad), t(bad, 1));

badIntervals = ConsolidateIntervals(sortrows([artefactInterval; noisyInterval]));
ok = ~isnan(values);
clean = [t(ok), values(ok)];