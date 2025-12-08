function varargout = PETH(samples, events, varargin)

% PETH - Compute a peri-event time histogram relative to synchronizing events
%
%  USAGE
%
%    [mat,t,mean] = PETH(samples,events,<options>)
%
%  INPUT
%
%    samples         either a list of timestamps (e.g. a spike train) or a
%                    (n,2) matrix of [timestamps values], in case of a continuous
%                    signal (e.g. reactivation strength, local field potential)
%                    
%    events          timestamps to synchronize on (e.g., brain
%                    stimulations)
%    <options>       optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties     Values
%    -------------------------------------------------------------------------
%     'durations'    durations before and after synchronizing events for each
%                    trial (in s) (default = [-1 1])
%     'nBins'        number of time bins around the events (default = 101)
%     'fast'         if 'off' (default), sort 'samples' and 'events' before
%                    operating, otherwise they are expected to be sorted
%                    (only for timestamps input)
%     'group'        indeces to group samples, a separate output will be computed
%                    per group (much faster than calling PETH once per group;
%                    only for timestamps input)
%     'mode'         whether the sample data is linear ('l') or circular ('c')
%                    (for example, in the case 'samples' is the phase of an
%                    oscillation; only for matrix input)
%     'show'         display the mean PETH (default = 'on' when no outputs are
%                    requested and 'off' otherwise)
%     'smooth'       standard deviation for Gaussian kernel (default = 1 bin)
%                    applied to the mean peri-event activity 'm' (note, no
%                    smoothing will be applied to the output 'matrix')
%     'title'        if the results are displayed ('show' = 'on'), specify a
%                    desired title (default is deduced by variable names)
%     <plot options> any other property (and all the properties that follow)
%                    will be passed down to "plot" (e.g. 'r', 'linewidth', etc);
%                    because all the following inputs are passed down to "plot",
%                    make sure you put these properties last
%    =========================================================================
%
%  OUTPUT
%
%    mat             matrix containing the counts of a point process (for 
%                    timestamp data) or the avarage activity (for a continous
%                    signal) around the synchronizing events. Each column
%                    corresponds to a particular delay around the event (delay
%                    value indicated in timeBins), and each row corresponds to
%                    a synchronizing event ('mat' is a cell array when 'group'
%                    is provided)
%    timeBins        time bin delay values corresponding columns of 'mat'
%    mean            average activity across all events (a cell array when 'group' 
%                    is provided)
%
%  EXAMPLE
%
%    % show mean spiking activity around the stimuli:
%    PETH(spikes(:,1),stimuli); 
%
%    % compute the mean lfp around delta wave peaks:
%    [matrix,timeBins,m] = PETH([lfp.timestamps double(lfp.data(:,1))],deltaWaves.peaks); 
%
%    % get the order of delta wave power:
%    [~,order] = sort(deltaWaves.peakNormedPower); 
%    plot the mean lfp around delta waves as ordered according to delta wave power
%    PlotColorMap(matrix(order,:),'x',timeBins);
%
%  SEE
%
%    See also Sync, SyncHist, SyncMap, PlotSync, PETHTransition.

% Copyright (C) 2018-2025 by Ralitsa Todorova & Michaël Zugaro, Pietro Bozzo
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% default values
duration = [-1 1];
nBins = 101;
fast = 'off';
group = [];
mode = 'l';
smooth = 1;
namestring = replace([inputname(1) ' synchronised to ' inputname(2)],'_','\_');
pictureoptions = {};
if nargout < 1
    show = 'on';
else
    show = 'off';
end

for i = 1 : 2 : length(varargin)
    switch lower(varargin{i})
        case 'durations'
            duration = varargin{i+1};
            if ~isvector(duration) || length(duration) ~= 2
                error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'duration'
            duration = varargin{i+1};
            if ~isvector(duration) || length(duration) ~= 2
                error('Incorrect value for property ''durations'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'nbins'
            nBins = varargin{i+1};
            if ~isvector(nBins) || length(nBins) ~= 1
                error('Incorrect value for property ''nBins'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'fast'
            fast = varargin{i+1};
            if ~isastring(fast,'on','off')
                error('Incorrect value for property ''fast'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'group'
            group = varargin{i+1};
            if ~isdvector(group)
                error('Incorrect value for property ''group'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'show'
            show = varargin{i+1};
            if ~isastring(show,'on','off')
                error('Incorrect value for property ''show'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'mode'
            mode = varargin{i+1};
            if ~isastring(mode,'l','c')
                error('Incorrect value for property ''mode'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'title'
            namestring = varargin{i+1};
            if ~isastring(namestring)
                error('Incorrect value for property ''title'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        case 'smooth'
            smooth = varargin{i+1};
            if ~isvector(smooth) || length(smooth) ~= 1
                error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help PETH">PETH</a>'' for details).');
            end
        otherwise
            pictureoptions = varargin(i:end); break
    end
end

if isempty(group)
    groups = 1;
else
    groups = unique(group);
end

if size(samples,2) == 2
    % samples is a signal
    t = linspace(duration(:,1),duration(2),nBins);
    mat_t = bsxfun(@plus,events,t);
    dt = diff(samples(:,1));
    samples(dt>median(dt)*2,2) = nan;       % To take care of gaps in the signal : interpolate values in the gaps to nans
    if strcmp(mode,'l')
        mat = interp1(samples(:,1),samples(:,2),mat_t);
        m = Smooth(mean(mat,'omitnan'),smooth);
    else % circular data
        unwrapped = unwrap(samples(:,2));
        mat = wrap(interp1(samples(:,1),unwrapped,mat_t));
        angles = mean(exp(1i*mat),'omitnan');
        smoothed = [Smooth(imag(angles(:)),smooth) Smooth(real(angles(:)),smooth)];
        m = atan2(smoothed(:,1),smoothed(:,2));
    end
else
    % samples is a point process
    % synchronize samples to events
    [sync,Ie,Is] = Sync(samples,events,'durations',duration,'fast',fast);
    t = linspace(duration(1),duration(2),nBins+1); % nBins+1 chosen to match previous behavior of Bins
    time_bin = t(2) - t(1);
    % compute matrix
    if isempty(group)
        mat = sync2mat(sync,Ie,t,size(events,1),nBins);
        if strcmpi(show,'on') || nargout > 2
            m = smoothdata(mean(mat),'gaussian',5*smooth) / time_bin; % factor 5 chosen to match previous behavior of Smooth
        end
    else
        % compute per group
        group = group(Is);
        mat = cell(numel(groups),1);
        m = cell(numel(groups),1);
        for g = 1 : numel(groups)
            mat{g} = sync2mat(sync(group==g),Ie(group==g),t,size(events,1),nBins);
            if strcmpi(show,'on') || nargout > 2
                m{g} = smoothdata(mean(mat{g}),'gaussian',5*smooth) / time_bin;
            end
        end
    end
    % adjust times and title
    t = (t(1:end-1) + t(2:end)) / 2;
    if strcmpi(show,'on')
        namestring = [namestring ', ' num2str(numel(Ie)) ' x ' num2str(numel(unique(Ie))) ' instances'];
    end
end

% plot
if strcmpi(show,'on')
    if isempty(group)
        plot(t,m,pictureoptions{:});
    else
        hold on
        cellfun(@(x) plot(t,x,pictureoptions{:}),m)
    end
    title(namestring);
end

if nargout > 0
    varargout{1} = mat;
    varargout{2} = t;
end
if nargout > 2
    varargout{3} = m;
end

end

% --- helper functions ---

function mat = sync2mat(sync,Ie,t,nEvents,nBins)

    mat = zeros(nEvents,nBins);
    if ~isempty(sync)
        s = discretize(sync,t);
        mat(:) = accumarray(sub2ind(size(mat),Ie,s),1,[numel(mat),1]); % Can maybe change with size(mat) !!
    end

end