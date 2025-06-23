function ratio = FiringRatio(spikes,opt)

%FiringRatio - Compute ratio of active units of a population over time.
%
% Compute instantaneous number of firing units over the total number of
% units in a population.
%
%  USAGE
%
%    ratio = FiringRatio(spikes,<options>)
%
%    spikes         spike train, each row is [spike time, unit id]
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default is 0.05)
%     'limits'      [start stop] in s (default is ~ first and last spike times)
%     'smooth'      gaussian kernel std in number of samples, default is 1
%                   (no smoothing)
%    =========================================================================
%
%  OUTPUT
%
%    ratio          each row is [t, istantaneuous ratio]
%
%  SEE
%
%    See also Frequency.

% Copyright (C) 2025 by Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

arguments
    spikes (:,2) {mustBeNumeric}
    opt.binSize (1,1) {mustBeNumeric,mustBePositive} = 0.05
    opt.limits (2,1) {mustBeNumeric} = [NaN,NaN]
    opt.smooth (1,1) {mustBeNumeric,mustBePositive} = 1
end

% set default value
default_limits = [spikes(1,1)-10*opt.binSize, spikes(end,1)+10*opt.binSize];
opt.limits(isnan(opt.limits)) = default_limits(isnan(opt.limits));

% keep spikes inside limits and shift them to 0
spike_times = spikes(:,1);
spike_times = spike_times(spike_times >= opt.limits(1) & spike_times <= opt.limits(2));
spike_times = spike_times - opt.limits(1) + opt.binSize/2;

% discretize time
bin_ind = ceil(spike_times / opt.binSize);
bin_ind(bin_ind==0) = 1;

% ratio
n_units = numel(unique(spikes(:,2)));
ratio = accumarray(bin_ind,spikes(:,2),[],@(x) numel(unique(x))) / n_units;
time = (opt.limits(1) : opt.binSize : opt.limits(2)+opt.binSize/2).';
ratio = [time, [ratio;zeros(numel(time)-numel(ratio),1)]];

% smoothing
ratio(:,2) = smoothdata(ratio(:,2),'gaussian',opt.smooth);