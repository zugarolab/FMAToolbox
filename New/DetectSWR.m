function ripples = DetectSWR(Channels, varargin)

%DetectSWR - Detect Sharp Wave Ripple in epochs of continuous data.
%
% DetectSWR     Detect Sharp Wave Ripple in epochs of continuous data. This
%               code is designed for linear probes, with at least one shank
%               that simultaneously records across the pyramidal cell layer
%               i.e. the superficial contacts record the ripple and at
%               least the deepest contact records the Sharp Wave (see the
%               .png files within the bz_DetectSWR for examples).
%                   |o  |
%                   |  o| <-- superficial channels record ripple
%                   |o  |
%                   |  o|
%                   |o  |
%                   |  o|
%                   |o  |
%                   |  o|  <-- deep channels record Sharp Wave
%                   \   /
%                    \ /
%                     .
%               Before running this reocord inspect your .lfp or .lfp file
%               to find that shank which best shows the sharp wave/ripple
%               profile across the shank e.g. Channels = 25:32.
%
% USE CASES
% 1) uses default parameters
% detected_swr = bz_DetectSWR( Filebase, Channels, Epochs )
%
% 2) supply parameters as a structure with fields an values
% detected_swr = bz_DetectSWR( ..., options )
%
% 3) supply parameters as Property/Value pairs
% detected_swr = bz_DetectSWR( ..., 'swBP', [2 50], 'WinSize', 300, etc. )
%   Filtering routines adapted from 'detect_hfos' by Eran Stark
%
% NOTE: -see MEMORY CHECK for how this routine handles potential memory
%       overflows associated with reading in large field potential files
%       a.k.a how I try to avoid running over into swap files
%
%       -If you would like to see how the detector performs against a hand
%       labeled .evt file you created: name the .evt file Filebase.rip.evt
%       with events 'start' and 'stop'.
%
%       -To improve detection, statistical thresholds are set using a local
%       window around each putative detection. This means global changes in
%       SWR and ripple magnitude aren't documented. Though, in exchange, this
%       routine is robust to non-stationarities in the data.
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com
% Convert to buzcode format: Andrea Navas-Olive, 2019
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUTS AND FREE PARAMETERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Channels:     an array of two channel IDs following neuroscope (base 0)
%               First RIPPLE channel, and second the SHARP WAVE channel.
%               A noise channel is the third channel (optional) - HLarsson
%                                                                   05/23
%               Detection in based on these order, please be careful
%               E.g. [1 11]
%
% All following inputs are optional:
%
% basepath:     full path where session is located (default pwd)
%               e.g. /mnt/Data/buddy140_060813_reo/buddy140_060813_reo
%
% Epochs:       a list of time stamps [start stop] in seconds demarcating
%               epochs in which to detect SWR. If this argument is empty,
%               the entire .lfp/.lfp file will be used.
%
% saveMat:      logical (default=false) to save in buzcode format
%
% forceDetect   true or false to force detection (avoid load previous
%               detection, default false)
%
% swBP:         [low high], passband for filtering sharp wave activity
%
% ripBP:        [low high] passband for filtering around ripple activity
%
% per_thresswD: a threshold placed upon the sharp wave difference magnitude
%               of the candidate SWR cluster determined via k-means.
% per_thresRip: a threshold placed upon ripple power based upon the non-SWR
%               cluster determined via k-means.
%
% WinSize:      window size in milliseconds for non-overlapping detections
%
% Ns_chk:       sets a window [-Ns_chk, Ns_chk] in seconds around a
%               candidate SWR detection for estimating local statistics.
%
% thresSDswD:   [low high], high is a threshold in standard deviations upon
%               the detected maximum sharp wave difference magnitude, based
%               upon the local distribution. low is the cutoff for the
%               feature around the detected peak for determining the
%               duration of the event, also based upon the local
%               distribution.
%
% thresSDrip:   [low high], high is a threshold in standard deviations upon
%               the detected maximum ripple magnitude, based
%               upon the local distribution. low is the cutoff for the
%               feature around the detected peak for determining the
%               duration of the event, also based upon the local
%               distribution.
% minIsi:       a threshold setting the minimum time between detections in
%               seconds.
%
% minDurSW:     a threshold setting the minimum duration of a localized
%               sharp-wave event.
%
% maxDurSW:     a threshold setting the maximum duration of a localized
%               sharp-wave event.
%
% minDurRP:     a threshold setting the minimum duration of a localized
%               ripple event associated with a sharp-wave.
%
% check:        boolean, a flag that triggers whether or not you wish to
%               view a figure with the sharpwave/ripple score and the default
%               k-means separation and refine it manually if needed (offers the
%               option to enter a manually curated cleaning/scoring step or
%               to enter debug mode) e.g. to ignore outliers. Note: this
%               option will require the user interaction.
%
% EVENTFILE:    boolean, a flag that triggers whether or not to write out
%               an event file corresponding to the detections (for
%               inspection in neuroscope).
%
% FIGS:         boolean, a flag that allows the user to output a figure
%               detailing the final clustering of the features and SWR
%
% TRAINING:     boolean, a flag for processing labeled training data
%               supplied by the user.
%
% DEBUG:        boolean, a flag that triggers the program to output an
%               additional set of figures and information.
%
% MAXGIGS:      the maximum number of GIGABYTES this
%
% noPrompts     true/false disable any user prompts (default: true)
%
% useSPW        Reject all events without a significant (and sufficiently long)
%               sharp-wave true, and determine the [start stop] boundaries of
%               the final events exclusively using the sharp wave (default=true)
%
%%%%%%%%%%%%%%%
%%% OUTPUTS %%%
%%%%%%%%%%%%%%%
% detected_swr: structure with fields
%
%        .timestamps : Nx2, Start and end timestamps of each SWR
%        .peaks : Nx1, trough nearest to peak power of ripple timestamp
%                 of each SWR
%        .peakNormedPower : Peak power values of each SWR
%        .detectorName : Name of this function
%        .stdev : Standard deviation used to threshold ripple power
%        .SwMax : Sharp-wave magnitude [zscore "local empirical percentile"]
%        .RipMax : Ripple magnitude [zscore "local empirical percentile"]
%        .detectorinfo : a structure with fields corresponding to all the
%               parameters used in the analysis
%
% OPTIONAL
% event file:   an event file formatted as filename.R01.evt that is
%               compatible with neuroscope for inspecting detections
%               [peak start stop]
%
% The analysis also writes out log file (.txt) containing the code used.
%
% DEPENDENCIES: all the major routines are contained within the associated
%               private folder. All should need to do is add the bz_DetectSWR
%               folder to your path. One of the training processes runs
%               linear discriminant analysis on the joint and marginal
%               spaces of the features. There is a statistics toolbox
%               dependency, R2013a, which I believe requires different
%               commands in later versions. I can update this if this
%               routine becomes widely used.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SET DEFAULT FREE PARAMETERS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addParameter(p, 'basepath', pwd, @isfolder);
addParameter(p, 'Epochs', [], @isnumeric);
addParameter(p, 'saveMat', true, @islogical);
addParameter(p, 'forceDetect', false, @islogical);
addParameter(p, 'swBP', [2, 50], @isnumeric);
addParameter(p, 'ripBP', [80, 250], @isnumeric);
addParameter(p, 'per_thresswD', 10, @isnumeric); % percentile
addParameter(p, 'per_thresRip', 50, @isnumeric); % percentile
addParameter(p, 'WinSize', 40, @isnumeric); % ms
addParameter(p, 'Ns_chk', 5, @isnumeric); % seconds
addParameter(p, 'thresSDswD', [0.5, 2.5], @isnumeric); % standard deviations
addParameter(p, 'thresSDrip', [0.5, 2.5], @isnumeric); % standard deviations
addParameter(p, 'minIsi', 0.05, @isnumeric); % seconds
addParameter(p, 'minDurSW', 0.02, @isnumeric); % seconds
addParameter(p, 'maxDurSW', 0.5, @isnumeric); % seconds
addParameter(p, 'minDurRP', 0.025, @isnumeric); % seconds
addParameter(p, 'check', true, @islogical);
addParameter(p, 'EVENTFILE', true, @islogical);
addParameter(p, 'FIGS', false, @islogical);
addParameter(p, 'TRAINING', false, @islogical);
addParameter(p, 'DEBUG', false, @islogical);
addParameter(p, 'MAXGIGS', 16, @isnumeric); % GigaBytes
addParameter(p, 'noPrompts', true, @islogical);
addParameter(p, 'useSPW', true, @islogical);
addParameter(p, 'useEEG', false, @islogical); % if you want the function to look for a .eeg file before looking for the .lfp

parse(p, varargin{:});

basepath = p.Results.basepath;
Epochs = p.Results.Epochs;
saveMat = p.Results.saveMat;
forceDetect = p.Results.forceDetect;
swBP = p.Results.swBP;
ripBP = p.Results.ripBP;
per_thresswD = p.Results.per_thresswD;
per_thresRip = p.Results.per_thresRip;
WinSize = p.Results.WinSize;
Ns_chk = p.Results.Ns_chk;
thresSDswD = p.Results.thresSDswD;
thresSDrip = p.Results.thresSDrip;
minIsi = p.Results.minIsi;
minDurSW = p.Results.minDurSW;
maxDurSW = p.Results.maxDurSW;
minDurRP = p.Results.minDurRP;
check = p.Results.check;
EVENTFILE = p.Results.EVENTFILE;
FIGS = p.Results.FIGS;
TRAINING = p.Results.TRAINING;
DEBUG = p.Results.DEBUG;
MAXGIGS = p.Results.MAXGIGS;
noPrompts = p.Results.noPrompts;
useSPW = p.Results.useSPW;
useEEG = p.Results.useEEG;

%%%%%%%%%%%%%%%%%%%%
%%% INPUT CHECKS %%%
%%%%%%%%%%%%%%%%%%%%
mfname = mfilename;
% 1) fileabase
% user might have input an absolute or relative path
Filebase = basenameFromBasepath(basepath);
Filebase = fullfile(basepath, Filebase);
[pathname, filename, extname] = fileparts(Filebase);
if isempty(pathname)
    pathname = pwd;
end

if check, useSPW = false; end
% check whether user provided a file extension
if isempty(extname)
    filechk = dir([Filebase, '.xml']);
    if isempty(filechk)
        warning(['%s: incompatible format for 1st argument (Filebase)\n', ...
            '\tIf user supplied a relative path, the .xml and all\n', ...
            '\trelevant files must be in the current directory.\n'], mfname);
    end
else
    filechk = dir(basepath);
    if isempty(filechk)
        warning(['%s: incompatible format for 1st argument (Filebase)\n', ...
            '\tIf user supplied a relative path, the .xml and all\n', ...
            '\trelevant files must be in the current directory.\n'], mfname);
    end
end

% Read in xml of file parameters
if exist([Filebase, '.xml'], 'file')
    fprintf(1, 'Loading XML file...\n');
    par = LoadXml([Filebase, '.xml']);
    % pull parameters from .xml file
    SR = par.lfpSampleRate; % lfp sampling rate
    nChan = par.nChannels; % number of channels in the recording
else
    fprintf(1, 'Loading session file...\n')
    load([Filebase, '.session.mat'], 'session')
    SR = session.extracellular.srLfp;
    nChan = session.extracellular.nChannels;
    % extract precision
    par.nBits = cellfun(@str2double, regexp(session.extracellular.precision, '\d+', 'match'));
end

% 2) Channels
if ~all(Channels+1 > 0 & Channels+1 <= nChan) % corrected by AFR
    error(['%s: incompatible input for 2nd argument (Channels)\n', ...
        '\tAccording to the .xml file, the user did not supply\n', ...
        '\ta valid list of channels.\n'], mfname);
end
Nchan = length(Channels);

% Access lfp looking for an .lfp file first and a .lfp file second
if ~isempty(dir([Filebase, '.lfp'])) && ~useEEG
    lfp_file = [Filebase, '.lfp'];
    lfp_info = dir(lfp_file);
elseif ~isempty(dir([Filebase, '.eeg']))
    lfp_file = [Filebase, '.eeg'];
    lfp_info = dir(lfp_file);
else
    error(['%s: Field potential file could not be found\n', ...
        '\tThere is neither a .lfp nor a .lfp file matching the\n', ...
        '\tbase name of the directory supplied by the user.\n'], mfname);
end

% 3) Epochs
if nargin < 3 || isempty(Epochs)
    % use the whole recording if no restrictions supplied by user
    fprintf(1, '%s: Using Default Epochs:\n', mfname);
    fprintf(1, '\tNo Epochs supplied by user. Using the whole recording.\n');
    Epochs = [1 / SR, lfp_info.bytes / (2 * nChan * SR)];
elseif ~all(Epochs(:, 1) < Epochs(:, 2)) || ...
        ~all(Epochs(:) > 0 & Epochs(:) < lfp_info.bytes/(2 * nChan * SR)) || ...
        size(Epochs, 2) > 2
    error(['%s: Epochs are not formatted properly:\n', ...
        '\tEpochs must be formatted as a list of [start stop] times in seconds.\n'], mfname);
end
Nepochs = size(Epochs, 1);


%%%%%%%%%%%%%%%%%%%%
%%% MEMORY CHECK %%%
%%%%%%%%%%%%%%%%%%%%
% This routine requires 4*NumberOfChannelsRequested + 2 GigaBytes.
% This estimate includes memory used to copy the inputs into functions (my
% kingdom for a pointer!).
if isfield(par, 'nBits')
    ByteConversion = 64 / par.nBits; % (double in Bytes)/(nBits -> Bytes)
else
    ByteConversion = 4;
end
BytesPerChan = (lfp_info.bytes / nChan) * ByteConversion;
GigsRequested = (4 * BytesPerChan * Nchan + 2 * BytesPerChan) / (1024^3);
if GigsRequested > MAXGIGS
    error(['%s: Too much memory requested\n', ...
        '\tThe user requested %4.1f GigaBytes and the maximum\n', ...
        '\tcurrently allowed is set by the variable MAXGIGS at %4.1f GigaBytes.\n', ...
        '\tIf your system can handle more memory, increase the value of MAXGIGS.\n'], mfname, GigsRequested, MAXGIGS);
end

%%%%%%%%%%%%%%%%%%%%

if TRAINING
    % Access hand labeled ripple events and convert to seconds
    fprintf(1, 'Loading Labeled Event File for Training...\n');

    training_rip = [Filebase, '.rip.evt'];
    if isempty(dir(training_rip))
        error(['%s: Training Event file could not be found\n', ...
            '\tThere is not an event file "Filebase.rip.evt" in \n', ...
            '\tthe directory supplied by the user.\n'], mfname);
    end
    ripfile = fopen(training_rip, 'r');

    % count the number of lines in the file and rewind
    count = 0;
    while ischar(fgetl(ripfile))
        count = count + 1;
    end
    frewind(ripfile);

    % Access ripple events and convert to seconds [start, stop]
    if mod(count, 2)
        fclose(ripfile);
        error('Check your event file: you should have start/stop pairs, but you have an odd number of lines');
    end

    % Parse event file to pull out hand labeled start/stop event pairs
    rippleTs = zeros(count/2, 2);
    count = 1;
    while 1
        tline = fgetl(ripfile);
        if ~ischar(tline), break, end

        if strfind(tline, 'start')

            ws = regexp(tline, '\s');
            rippleTs(round(count/2), 1) = str2double(tline(1:ws-1)) / 1000;

        elseif strfind(tline, 'stop')

            ws = regexp(tline, '\s');
            rippleTs(round(count/2), 2) = str2double(tline(1:ws-1)) / 1000;
        end
        count = count + 1;
    end
    fclose(ripfile);

    % For each ripple, I will select a segment of lfp data of the same length,
    % that does not include a ripple. I will use these paired samples for
    % comparison.
    nonrippleTs = zeros(size(rippleTs));
    counter = 0;

    % for each trial, find the labeled ripples contained
    fprintf(1, 'Parsing non-ripple paired samples...\n');
    for ii = 1:Nepochs

        ripples = find(rippleTs(:, 1) > Epochs(ii, 1) & rippleTs(:, 1) < Epochs(ii, 2));
        if isempty(ripples)
            continue
        end

        % for each ripple within this trial, select a corresponding nonripple
        % segment, checking to make sure it doesn't overlap with any ripples
        % within that trial.
        Trip = rippleTs(ripples, :);
        Nrip = length(ripples);
        nonripTmp = zeros(Nrip, 2);
        nonripTmp(:, 1) = (Epochs(ii, 2) - Epochs(ii, 1)) * rand(Nrip, 1) + Epochs(ii, 1);
        nonripTmp(:, 2) = nonripTmp(:, 1) + diff(Trip, [], 2);

        % Check nonripple segments for overlap with ripple segments and each
        % other
        flags = true(Nrip, 1);
        inds = 1:Nrip;
        while any(flags)
            for jj = inds
                if flags(jj)

                    % Concatenate all nonripples and all ripple times,
                    % excluding the one under inspection.
                    if Nrip > 1
                        temp = nonripTmp;
                        temp(jj, :) = []; % exclude current nonripple
                        Temp = zeros(2*Nrip-1, 2);
                        Temp(1:Nrip-1, :) = temp;
                        Temp(Nrip:end, :) = Trip;
                    else
                        Temp = Trip;
                    end

                    % Check for overlap between current "nonripple" and ripples
                    % and other "nonripples"
                    if isempty(find(nonripTmp(jj, 1) > Temp(:, 1) & nonripTmp(jj, 1) < Temp(:, 2), 1)) && ...
                            isempty(find(nonripTmp(jj, 2) > Temp(:, 1) & nonripTmp(jj, 2) < Temp(:, 2), 1))

                        % log nonripple
                        nonrippleTs(counter+inds(jj), :) = nonripTmp(jj, :);
                        flags(jj) = 0;

                    else

                        % Resample
                        nonripTmp(jj, 1) = (Epochs(ii, 2) - Epochs(ii, 1)) * rand + Epochs(ii, 1);
                        nonripTmp(jj, 2) = nonripTmp(jj, 1) + diff(Trip(jj, :));

                        % set all flags to true to double check that this
                        % resample didn't result in overlap
                        flags = true(Nrip, 1);
                    end
                end
            end
        end
        % increment counter
        counter = counter + Nrip;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ESTIMATE LFP FEATURES: START %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% As a consession to memory usage, I interleave loading of field potential
% variables, and extracting segments in the case of training data.

% load lfp file
fprintf(1, 'Loading LFP file...\n');
lfp        = LoadBinary(lfp_file,'frequency',SR,'nChannels',nChan, ...
    'channels',Channels+1);

lfp = int16(lfp);

% filter for sharp wave using Eran's routines and compute SW diff feature
hSw1 = makegausslpfir(swBP(1), SR, 6);
hSw2 = makegausslpfir(swBP(2), SR, 6);

lfpLow = firfilt(lfp, hSw2); % lowpass filter
eegLo = firfilt(lfpLow, hSw1); % highpass filter
% lfpLow     = lfpLow -.lfpLo;            % difference of Gaussians
% clear .lfpLo                          %DL: not sure what this '.' was doing here...
lfpLow = lfpLow - eegLo; % difference of Gaussians
clear eegLo


swDiff = lfpLow(:, 1) - lfpLow(:, 2); % this allows us to use a noise channel

if TRAINING
    % Extract ripple and paired non-ripple segments from lfp
    fprintf(1, 'Training: Extracting Paired Segments of lfp, lfpLow, and swDiff...\n');

    % convert timestamps to samples
    rippleSamp = floor(rippleTs*SR);
    nonrippleSamp = floor(nonrippleTs*SR);

    % cut out ripple and paired non-ripple segments
    pairs0 = cell(counter, 2);
    pairs1 = cell(counter, 2);
    pairs2 = cell(counter, 2);

    for ii = 1:counter
        pairs0{ii, 1} = lfp(rippleSamp(ii, 1):rippleSamp(ii, 2), :);
        pairs0{ii, 2} = lfp(nonrippleSamp(ii, 1):nonrippleSamp(ii, 2), :);

        pairs1{ii, 1} = lfpLow(rippleSamp(ii, 1):rippleSamp(ii, 2), :);
        pairs1{ii, 2} = lfpLow(nonrippleSamp(ii, 1):nonrippleSamp(ii, 2), :);

        pairs2{ii, 1} = swDiff(rippleSamp(ii, 1):rippleSamp(ii, 2), 1);
        pairs2{ii, 2} = swDiff(nonrippleSamp(ii, 1):nonrippleSamp(ii, 2), 1);

    end
end
% clear lfpLow % this can be useful if the user wants to manually remove ripples with positive sharp waves while in mode "check"

% I made 2 changes to Eran's defaults for ripple processing:
% 1) I subtract the across channel mean from each channel prior to
% filtering
% 2) After subtracting, filtering, envelope, powerwin, and take the maximum
% across all channels at each time point as the ripple power.

% filter and process ripples, using Eran's defaults (diff of Gaussians)
hRip1 = makegausslpfir(ripBP(1), SR, 6);
hRip2 = makegausslpfir(ripBP(2), SR, 6);
if Nchan == 2
    % Use the difference between the two channels for computing ripple power
    % subtract mean before filtering
    lfpM = lfp - int16(repmat(mean(lfp, 2), [1, Nchan]));
    %clear lfp

    rip = firfilt(lfpM, hRip2); % highpass filter
    clear lfpM
    eegLo = firfilt(rip, hRip1); % lowpass filter
    % rip        = rip -.lfpLo;               % difference of gaussians
    % clear.lfpLo
    rip = rip - eegLo; % difference of gaussians
    clear eegLo

    ripWindow = pi / mean(ripBP);
    powerWin = makegausslpfir(1/ripWindow, SR, 6);

    rip = abs(rip);
    ripPower0 = firfilt(rip, powerWin);
    clear rip

    ripPower0 = max(ripPower0, [], 2);


else % if a noise channel is provided

    % ripple channel
    rip = firfilt(lfp(:, 1), hRip2); % highpass filter
    eegLo = firfilt(rip, hRip1); % lowpass filter
    rip = rip - eegLo; % difference of gaussians
    ripWindow = pi / mean(ripBP);
    powerWin = makegausslpfir(1/ripWindow, SR, 6);
    rip = abs(rip);
    signal = firfilt(rip, powerWin);

    rip = firfilt(lfp(:, 3), hRip2); % highpass filter
    eegLo = firfilt(rip, hRip1); % lowpass filter
    rip = rip - eegLo; % difference of gaussians
    ripWindow = pi / mean(ripBP);
    powerWin = makegausslpfir(1/ripWindow, SR, 6);
    rip = abs(rip);
    noise = firfilt(rip, powerWin);

    ripPower0 = signal - noise;
end


if TRAINING
    % Extract ripple and paired non-ripple segments from lfp
    fprintf(1, 'Training: Extracting Paired Segments of ripple Power...\n');

    % cut out ripple and paired non-ripple segments
    pairs3 = cell(counter, 2);

    for ii = 1:counter

        pairs3{ii, 1} = ripPower0(rippleSamp(ii, 1):rippleSamp(ii, 2), 1);
        pairs3{ii, 2} = ripPower0(nonrippleSamp(ii, 1):nonrippleSamp(ii, 2), 1);

    end

    if DEBUG
        % Inspect paired segements
        figure('Position', get(0, 'ScreenSize'))
        % Examine pairs, including low pass filtered data
        for ii = 1:counter
            set(gcf, 'Name', sprintf('SWR Pair %d', ii))
            for jj = 1:Nchan

                h1 = subplot(Nchan+2, 2, 2*jj-1); ...
                    plot(pairs0{ii, 1}(:, jj)), hold on
                plot(pairs1{ii, 1}(:, jj), 'r'), hold off
                ylims1 = get(h1, 'ylim');

                h2 = subplot(Nchan+2, 2, 2*jj); ...
                    plot(pairs0{ii, 2}(:, jj)), hold on
                plot(pairs1{ii, 2}(:, jj), 'r'), hold off
                ylims2 = get(h2, 'ylim');

                if jj == 1
                    set(get(h1, 'title'), 'string', 'Labeled SWR')
                    set(get(h2, 'title'), 'string', 'Paired Non-SWR sample')
                end
                set(h1, 'ylim', [min([ylims1(1), ylims2(1)]), max([ylims1(2), ylims2(2)])])
                set(h2, 'ylim', [min([ylims1(1), ylims2(1)]), max([ylims1(2), ylims2(2)])])
            end

            h3 = subplot(Nchan+2, 2, 2*Nchan+1);
            plot(pairs2{ii, 1}, 'b'), hold off
            ylims3 = get(h3, 'ylim');
            ylabel('SW Diff')
            h4 = subplot(Nchan+2, 2, 2*Nchan+2);
            plot(pairs2{ii, 2}, 'r'), hold off
            ylims4 = get(h4, 'ylim');
            h5 = subplot(Nchan+2, 2, 2*(Nchan + 1)+1);
            plot(pairs3{ii, 1}, 'b'), hold off
            ylims5 = get(h5, 'ylim');
            ylabel('Ripple Power')
            h6 = subplot(Nchan+2, 2, 2*(Nchan + 1)+2);
            plot(pairs3{ii, 2}, 'r'), hold off
            ylims6 = get(h6, 'ylim');

            set(h3, 'ylim', [min([ylims3(1), ylims4(1)]), max([ylims3(2), ylims4(2)])])
            set(h4, 'ylim', [min([ylims3(1), ylims4(1)]), max([ylims3(2), ylims4(2)])])
            set(h5, 'ylim', [min([ylims5(1), ylims6(1)]), max([ylims5(2), ylims6(2)])])
            set(h6, 'ylim', [min([ylims5(1), ylims6(1)]), max([ylims5(2), ylims6(2)])])

            drawnow
            pause(.1)
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ESTIMATE LFP FEATURES: END %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TRAINING

    fprintf(1, 'Estimating Performance on Training Data...\n');

    % Assess how informative these features are for detection
    % as a first pass, I will take the max of each feature with each paired
    % sample window
    swDiffMax = zeros(counter, 2);
    ripPowerMax = zeros(counter, 2);

    for ii = 1:counter
        swDiffMax(ii, 1) = max(pairs2{ii, 1});
        swDiffMax(ii, 2) = max(pairs2{ii, 2});

        ripPowerMax(ii, 1) = max(pairs3{ii, 1});
        ripPowerMax(ii, 2) = max(pairs3{ii, 2});
    end

    % Calculate marginal histograms
    [n1a, ctr1a] = hist(swDiffMax(:, 1), 40);
    [n1b, ctr1b] = hist(swDiffMax(:, 2), 40);
    [n2a, ctr2a] = hist(ripPowerMax(:, 1), 40);
    [n2b, ctr2b] = hist(ripPowerMax(:, 2), 40);

    % run a basic k-means clustering and evaluate performance
    [idx, Ctrs] = kmeans([swDiffMax(:), ripPowerMax(:)], 2);

    labels = [ones(counter, 1); 2 * ones(counter, 1)];
    scores = logical(idx == labels);

    % let's run a linear discriminant analysis on the joint and marginal spaces
    % (requires the statistics toolbox, R2013a, I believe the commands change
    % in later versions)
    FeaturesAll = [swDiffMax(:), ripPowerMax(:)];
    cls1 = ClassificationDiscriminant.fit(FeaturesAll, labels); % joint
    cls2 = ClassificationDiscriminant.fit(FeaturesAll(:, 1), labels); % swDiff
    cls3 = ClassificationDiscriminant.fit(FeaturesAll(:, 2), labels); % rip Power
    % access performance on training data
    perf0 = 100 - ((resubLoss(cls1) * cls1.NObservations) / cls1.NObservations) * 100;
    perf1 = 100 - ((resubLoss(cls2) * cls2.NObservations) / cls2.NObservations) * 100;
    perf2 = 100 - ((resubLoss(cls3) * cls3.NObservations) / cls3.NObservations) * 100;

    % For the unsupervised case, I'll set thresholds for ripple power
    % and sharp wave difference magnitude according to thresholds upon the
    % empirical distributions. This requires I map the probabilities of one
    % cluster to that of the other.
    x = ripPowerMax(:, 1);
    [Fi, xi] = ecdf(x);
    n = length(x);
    xj = xi(2:end);
    Fj = (Fi(1:end-1) + Fi(2:end)) / 2;
    xj = [xj(1) - Fj(1) * (xj(2) - xj(1)) / ((Fj(2) - Fj(1))); ...
        xj; ...
        xj(n) + (1 - Fj(n)) * ((xj(n) - xj(n-1)) / (Fj(n) - Fj(n-1)))];
    Fj = [0; Fj; 1];

    if DEBUG
        figure('Name', 'Empirical Cumulative Probability Distribution')
        stairs(xi, Fi, 'r');
        hold on
        plot(xj, Fj, 'b-');
        hold off
    end

    % Map the Nth percentile in the non-SWR cluster to the Pth in the SWR
    % cluster.
    per_thres = 50;
    per_corr = 100 * interp1(xj, Fj, prctile(ripPowerMax(:, 2), per_thres));

    % Main Figure in the case of hand labeled training data
    % close all
    figure('Color', 'w'),

    string = sprintf(['Sharp Wave Difference vs. Ripple Power\n', ...
        'Mean subtraction before filter\n', ...
        'Max at each time point is taken for ripple power']);

    % Joint space of SWR difference and ripple power
    h1 = subplot(2, 2, 2);
    plot(swDiffMax(:, 1), ripPowerMax(:, 1), '.b'), hold on
    plot(swDiffMax(:, 2), ripPowerMax(:, 2), '.r'), hold on
    axis([min(swDiffMax(:)), max(swDiffMax(:)), 0, max(ripPowerMax(:))])
    % LDA coefficients
    K1 = cls1.Coeffs(1, 2).Const;
    L1 = cls1.Coeffs(1, 2).Linear;
    % LDA function
    f1 = @(x1, x2) K1 + L1(1) * x1 + L1(2) * x2;
    % Percentile threshold
    line(get(h1, 'xlim'), repmat(prctile(ripPowerMax(:, 1), per_corr), [1, 2]), ...
        'linestyle', '--', 'color', 'k')
    l0 = ezplot(h1, f1, [get(h1, 'xlim'), get(h1, 'ylim')]);
    set(l0, 'Color', 'g', 'linewidth', 2)
    title(h1, string, 'fontsize', 12)
    legend('SWR', '~SWR', 'Location', 'NorthEastOutside')
    xlabel('Sharp Wave Diff', 'fontsize', 12)
    ylabel('Ripple Power', 'fontsize', 12)
    set(h1, 'fontsize', 12)
    % Marginal Space of Sharp Wave Difference
    h2 = subplot(2, 2, 4);
    bar(ctr1a, -n1a, 1, 'b'), hold on, bar(ctr1b, -n1b, 1, 'r'), ...
        axis([min(swDiffMax(:)), max(swDiffMax(:)), -max([n1a, n1b]) * 1.1, 0]),
    % LDA coefficients
    K2 = cls2.Coeffs(1, 2).Const;
    L2 = cls2.Coeffs(1, 2).Linear;
    % LDA function
    f2 = @(x1, x2) K2 + L2(1) * x1;
    l1 = ezplot(h2, f2, [get(h2, 'xlim'), get(h2, 'ylim')]);
    set(l1, 'Color', 'g', 'linewidth', 2)
    title('')
    xlabel('')
    ylabel('')
    axis('off'),
    % Marginal Space of Ripple Power
    h3 = subplot(2, 2, 1);
    barh(ctr2a, -n2a, 1, 'b'), hold on, barh(ctr2b, -n2b, 1, 'r'), ...
        axis([-max([n2a, n2b]) * 1.1, 0, min(ripPowerMax(:)), max(ripPowerMax(:))])
    % LDA coefficients
    K3 = cls3.Coeffs(1, 2).Const;
    L3 = cls3.Coeffs(1, 2).Linear;
    % LDA function
    f3 = @(x2, x1) K3 + L3(1) * x1; % note the flip in x and y
    l2 = ezplot(h3, f3, [get(h3, 'xlim'), get(h3, 'ylim')]);
    set(l2, 'Color', 'g', 'linewidth', 2)
    title('')
    xlabel('')
    ylabel('')
    axis('off'),
    set(h1, 'Position', [0.25, 0.35, 0.55, 0.55]);
    set(h2, 'Position', [.25, .1, .55, .15]);
    set(h3, 'Position', [.05, .35, .15, .55]);

    string = sprintf(['K-means Clustering:\n', ...
        '    Center1: X  = %7.2f, Y = %7.2f\n', ...
        '    Center2: X  = %7.2f, Y = %7.2f\n', ...
        '    Performance = %4.1f percent\n\n', ...
        'LDA, joint space:\n', ...
        '    Performance = %4.1f percent\n\n', ...
        'LDA, 1D SW Diff:\n', ...
        '    Performance = %4.1f percent\n\n', ...
        'LDA, 1D Ripple Power:\n', ...
        '    Performance = %4.1f percent\n\n', ...
        '%dth percentile of ~SWR corresponds to\n%3.1fth percentile of SWR'], ...
        Ctrs(1, 1), Ctrs(1, 2), Ctrs(2, 1), Ctrs(2, 2), ...
        (sum(scores) / numel(scores))*100, perf0, perf1, perf2, per_thres, per_corr);

    annotation('textbox', [0.81, 0.70, .1, .1], 'String', string, ...
        'edgecolor', 'none', 'fontsize', 12)

end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DETECT SWR EVENTS %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% I take the max of the sharp waves in non-overlapping windows. Once a
% maximum in the sharp wave is detected, the maximum ripple power around
% that detected is logged.

% Specifically, maximum detection is with reference to the Sharp Wave. If a
% 3 point maximum is detected for the sharp wave feature, then a maximum
% for the ripple power is also logged for that detection within a window
% around sharp wave detection.

% This first pass is permissive: it is just detecting maxima without
% regards for magnitude

WinSize = floor((WinSize * SR)/1000); %#ok<NODEF> % window in samples
HalfWinSize = floor(WinSize/2); % for logging ripple power around SW

% for each trial, block features into windows
% first pass, just count blocks for pre-allocating memory
Nblocks = 0;
for ii = 1:Nepochs
    epoch_inds = floor(Epochs(ii, :)*SR);
    blocks = epoch_inds(1):WinSize:epoch_inds(2);
    Nblocks = Nblocks + numel(blocks) - 1;
end

% Second pass populates values
swDiffAll = zeros(Nblocks, 1);
ripPowerAll = zeros(Nblocks, 1);

% In samples at SR
featureTs = zeros(Nblocks, 1);

MaxCounter = 0;
for ep_i = 1:Nepochs
    epoch_inds = floor(Epochs(ep_i, :)*SR);
    blocks = epoch_inds(1):WinSize:epoch_inds(2);

    for jj = 1:length(blocks) - 1
        %Find the maximum sharp-wave power in the block
        [tempswDiff, ind0] = max(swDiff(blocks(jj):blocks(jj+1)-1, 1));

        % Transform from sample within block to sample within record
        maxSample = blocks(jj) + ind0 - 1;

        % If the maximum is detected on the edges of the window, run a 3
        % point max check over the boundaries to resolve edge effect
        if (ind0 == 1 || ind0 == WinSize)

            % Check boundaries and then check for 3 point max
            if (maxSample == 1) || (maxSample == size(swDiff, 1))
                % Don't use the boundaries of the recording
                continue
            else
                % Check for 3 point maximum
                [~, maxchk] = max([swDiff(maxSample-1, 1), ...
                    swDiff(maxSample, 1), ...
                    swDiff(maxSample+1, 1)]);
                if abs(2-maxchk) > 0
                    % This is not a maximum
                    continue
                end
            end
        end
        % Log features
        MaxCounter = MaxCounter + 1;
        featureTs(MaxCounter) = maxSample;
        swDiffAll(MaxCounter) = tempswDiff;
        %ripPowerAll(MaxCounter) = max(ripPower0(maxSample-HalfWinSize:maxSample+HalfWinSize,1));
        %added by aza to detect very long
        if (maxSample - HalfWinSize) < 1
            ripPowerAll(MaxCounter) = max(ripPower0(1:maxSample+HalfWinSize, 1));
        else
            %aza
            try
                ripPowerAll(MaxCounter) = max(ripPower0(max([maxSample - HalfWinSize, 1]):min([maxSample + HalfWinSize, length(ripPower0)]), 1));
            catch
                ripPowerAll(MaxCounter) = max(ripPower0(maxSample-10:maxSample+10, 1));
            end
        end
    end
    featureTs = featureTs(1:MaxCounter);
    swDiffAll = swDiffAll(1:MaxCounter);
    ripPowerAll = ripPowerAll(1:MaxCounter);

    % run a basic k-means clustering to partition SWR from ~SWR
    [idx, Ctrs] = kmeans([swDiffAll(:), ripPowerAll(:)], 2);

    % we know there will be far fewer SWR than not, so let's assign idx1 to
    % the cluster with the fewest members i.e. SWR
    if sum(idx == 1) > sum(idx == 2)
        idx1 = logical(idx == 1);
        idx2 = logical(idx == 2);
        % reverse labels
        temp = idx1;
        idx1 = idx2;
        idx2 = temp;
        clear temp
    else
        idx1 = logical(idx == 1);
        idx2 = logical(idx == 2);
    end

    if check
        clf
        plot(swDiffAll(idx2), ripPowerAll(idx2), 'k.', 'markersize', 1); hold on
        plot(swDiffAll(idx1), ripPowerAll(idx1), 'r.', 'markersize', 1);
        legend('non-ripples', 'ripples');
        %         disp(['Please review the figure and change the idx1 (ripples) ' newline 'and idx2 (calm non-ripple periods) groups.' newline ...
        %             'Note that noisy periods need not participate in either group.']);
        choice = str2double(input('Would you like to perform the manual cleaning step? (1=yes; 0=no, keep groups as they are; -1=enter debug mode): ', 's'));
        if choice == 1
            DetectSWR_check_helper % this is a script in the "private" folder
        elseif choice == -1
            keyboard
        end
        idx2(swDiffAll > mean(swDiffAll(idx1))) = 0;
        idx2(ripPowerAll > mean(ripPowerAll(idx1))) = 0;
    end

    % For the unsupervised case, I'll set thresholds for ripple power
    % and sharp wave difference magnitude according to thresholds upon the
    % empirical distributions.

    if TRAINING
        fprintf(1, 'Estimating Precision and Recall surfaces for features...');
        % Construct Precision/Recall surfaces
        thresL_swDs = 0:5:100;
        thresL_Rips = 0:5:100;

        % Calculate precision/recall values if labeled data is supplied
        % precision: true_positive / (true_positive + false_positive)
        % recall:    true_positive / (true_positive + false_negative)
        Precision = zeros([length(thresL_swDs), length(thresL_Rips)]);
        Recall = zeros([length(thresL_swDs), length(thresL_Rips)]);

        for jj = 1:length(thresL_swDs)
            for kk = 1:length(thresL_Rips)

                per_thresswDt = thresL_swDs(jj);
                per_thresRipt = thresL_Rips(kk);

                % Sharp Wave Difference: set a lower cutoff in terms of percentiles
                thresL_swD = prctile(swDiffAll(idx1), per_thresswDt);

                % Map the Nth percentile in the non-SWR cluster to the Pth in the SWR
                % cluster.
                thresL_Rip = prctile(ripPowerAll(idx2, 1), per_thresRipt);

                % retain SWR candidates according to clustering and thresholds
                SWRind = idx1 & swDiffAll > thresL_swD & ripPowerAll > thresL_Rip;

                % convert featureTs to seconds
                featureTsS = featureTs * 1 / SR;
                SWRs = featureTsS(SWRind);

                TP = 0; % true positive
                FP = 0; % false positive

                hits = false(size(rippleTs, 1), 1);
                for ii = 1:length(SWRs)

                    chk = SWRs(ii) > rippleTs(:, 1) & SWRs(ii) < rippleTs(:, 2);

                    if sum(chk) < 1
                        FP = FP + 1;
                    else
                        % exclude duplicates
                        if any(hits & chk)
                            continue
                        end
                        TP = TP + 1;
                        hits = hits | chk;
                    end

                end

                FN = sum(~hits); % false negatives

                Precision(jj, kk) = TP / (TP + FP);
                Recall(jj, kk) = TP / (TP + FN);

            end
        end

        figure,
        subplot(121), imagesc(Precision),
        xlabel('thresL rip'), ylabel('thresL swD'), title('Precision'),
        set(gca, 'ytick', 1:2:21)
        set(gca, 'yticklabel', 0:10:100)
        set(gca, 'xtick', 1:2:21)
        set(gca, 'xticklabel', 0:10:100)
        colorbar
        axis square

        subplot(122), imagesc(Recall),
        xlabel('thresL rip'), ylabel('thresL swD'), title('Recall'),
        set(gca, 'ytick', 1:2:21)
        set(gca, 'yticklabel', 0:10:100)
        set(gca, 'xtick', 1:2:21)
        set(gca, 'xticklabel', 0:10:100)
        colorbar
        axis square
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LOCALIZING EVENTS %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Localize events and winnow detections using the following criteria:
    fprintf('Localizing and Screening Detections ...\n');
    % 1) detected sharp wave difference magnitude must be > thresSDswD(2)
    % 2) detected ripple power must be > thresSDrip(2)
    % 3) detections can't be within a minimum inter-event threshold
    % 4) both the detected sharp-wave and ripple power increase mubt persistent
    % for a minimum duration.
    minDurSWs = floor(minDurSW*SR); % convert to samples
    minDurRPs = floor(minDurRP*SR);
    % 5) The sharp-wave can not exceed a maximum duration
    maxDurSWs = floor(maxDurSW*SR);

    % Sharp Wave Difference: set a lower cutoff in terms of percentiles
    thresL_swD = prctile(swDiffAll(idx1), per_thresswD);

    % Set the value in the the non-SWR cluster at the Nth percentile as the
    % lower threshold for ripple detection
    thresL_Rip = prctile(ripPowerAll(idx2, 1), per_thresRip);

    % retain SWR candidates according to clustering and thresholds
    if check % if the points were clustered manually, there is no need to assign further thresholds on them
        SWRind = idx1;
    else
        SWRind = idx1 & swDiffAll > thresL_swD & ripPowerAll > thresL_Rip;
    end
    SWRid = find(SWRind);

    % Access time-stamps of SWR features
    SWRs = featureTs(SWRind);
    Nswr = length(SWRs);

    % Access the subset of detections in samples to check for localization
    bound = Ns_chk * SR; % convert to samples

    % Exclude candidate SWRs that don't have a +/- bounds around them. Of
    % course, if you make bounds huge, or have a short recording, this will be
    % a problem
    exclusions = (SWRs - bound < 1) | (SWRs + bound > size(swDiff, 1));
    SWRs(exclusions) = [];
    SWRind(exclusions) = [];
    SWRid(exclusions) = [];

    fprintf(1, 'Number of candidate SWRs excluded due to boundary constraints: %d\n', Nswr-length(SWRs));
    Nswr = length(SWRs);
    timeline = -bound:bound;
    SWR_diff = [diff([0; SWRs * (1 / SR)])]; % in seconds

    % localize events
    SWR_valid.Ts = zeros(size(SWRs, 1), 3); % peak, start, stop in samples
    SWR_valid.SwMax = zeros(size(SWRs, 1), 2); % z-score & percentile
    SWR_valid.RipMax = zeros(size(SWRs, 1), 2); % z-score & percentile

    valid = 0;
    SWRindF = SWRind;

    if DEBUG
        figure
    end

    for ii = 1:Nswr

        if DEBUG
            set(gcf, 'Name', sprintf('Detection # %d out of %d', ii, Nswr))
        end

        % threshold feature distributions to localize about detection

        % 1) detected sharp wave difference magnitude must be > thresSD locally
        swDiff_chk = swDiff(SWRs(ii)-bound:SWRs(ii)+bound);
        MswDiff = median(swDiff_chk);
        SDswDiff = std(swDiff_chk);

        if swDiff_chk(bound+1) < MswDiff + thresSDswD(2) * SDswDiff && useSPW
            if DEBUG
                disp('Sharp Wave Difference magnitude too small locally')
            end
            SWRindF(SWRid(ii)) = false;
            continue
        end

        % 2) detected ripple power must be > thresSD locally
        indices = (SWRs(ii) - bound:SWRs(ii) + bound)';
        ripPow_chk = [indices, ripPower0(SWRs(ii)-bound:SWRs(ii)+bound)];
        MripPow = median(ripPow_chk(:, 2));
        SDripPow = std(ripPow_chk(:, 2));

        if ripPowerAll(SWRid(ii)) < MripPow + thresSDrip(2) * SDripPow
            if DEBUG
                disp('Ripple Power too low locally')
            end
            SWRindF(SWRid(ii)) = false;
            continue
        end

        % 3) detections can't be within a minimum inter-event interval thresh
        if ii < Nswr
            if SWR_diff(ii) < minIsi
                if DEBUG
                    disp('IEI is too small')
                end
                SWRindF(SWRid(ii)) = false;
                continue
            end
        end

        % Estimate duration of detected Sharp-Wave event
        startSWR = find(swDiff_chk(1:bound+1) < MswDiff+thresSDswD(1)*SDswDiff, 1, 'last');
        stopSWR = find(swDiff_chk(bound+1:end) < MswDiff+thresSDswD(1)*SDswDiff, 1, 'first');
        durSWR = (bound + 1 - startSWR) + stopSWR;

        % Estimate duration of detected ripple power increase event
        % find index of maximum ripple power increase around Sharp-Wave
        potentialPeak = ripPow_chk(bound+1-HalfWinSize:bound+1+HalfWinSize, :);
        [~, maxRpInd] = max(potentialPeak(:, 2));
        potentialStart = ripPow_chk(1:bound-HalfWinSize+maxRpInd, :);
        startRP = potentialStart(find(potentialStart(:, 2) < MripPow+thresSDrip(1)*SDripPow, 1, 'last'), 1);
        potentialStop = ripPow_chk(bound-HalfWinSize+maxRpInd:end, :);
        stopRP = potentialStop(find(potentialStop(:, 2) < MripPow+thresSDrip(1)*SDripPow, 1, 'first'), 1);
        durRP = stopRP - startRP;

        % 4) events associated with detections must exceed a minimum duration
        if durRP < minDurRPs && durSWR < minDurSWs
            if DEBUG
                disp('Event Duration is too short.')
            end
            SWRindF(SWRid(ii)) = false;
            continue
        end

        % 5) sharp-wave must not exceed a maximum duration
        if durSWR > maxDurSWs && useSPW
            if DEBUG
                disp('Sharp-Wave Event Duration is too long.')
            end
            SWRindF(SWRid(ii)) = false;
            continue
        end

        % We have a valid detection
        valid = valid + 1;
        if useSPW, % let event boundaries be entirely determined by the sharp wave boundaries
            SWR_valid.Ts(valid, :) = [SWRs(ii), (SWRs(ii) - (bound + 1 - startSWR)), (SWRs(ii) + stopSWR - 1)];
        else % event boundaries are determined by the union of the sharp wave and ripple intervals
            SWR_valid.Ts(valid, :) = [SWRs(ii), min([SWRs(ii) - (bound + 1 - startSWR), startRP]), max([SWRs(ii) + stopSWR - 1, stopRP])];
        end

        swMax = swDiffAll(SWRid(ii));
        ripMax = ripPowerAll(SWRid(ii));

        SWR_valid.SwMax(valid, :) = ...
            [(swMax - MswDiff) / SDswDiff, sum(swDiff_chk < swMax) / (2 * bound + 1)];

        SWR_valid.RipMax(valid, :) = ...
            [(ripMax - MripPow) / SDripPow, sum(ripPow_chk(:, 2) < ripMax) / (2 * bound + 1)];

        if DEBUG
            h1 = subplot(2, 3, [1, 2]);
            plot(timeline, swDiff_chk), axis tight, hold on
            plot([timeline(startSWR), timeline(bound+stopSWR)], ...
                [swDiff_chk(startSWR), swDiff_chk(bound+stopSWR)], 'og')
            line([0, 0], get(gca, 'ylim'), ...
                'linestyle', '--', 'color', 'k')
            line(get(gca, 'xlim'), repmat(MswDiff, [1, 2]), ...
                'linestyle', '--', 'color', 'k')
            line(get(gca, 'xlim'), repmat(MswDiff+thresSDswD(2)*SDswDiff, [1, 2]), ...
                'linestyle', '--', 'color', 'r')
            line(get(gca, 'xlim'), repmat(MswDiff-thresSDswD(2)*SDswDiff, [1, 2]), ...
                'linestyle', '--', 'color', 'r')
            xlabel('Samples'),
            ylabel('A.U.'),
            title(sprintf(['Sharp Wave Difference Magnitude\n', ...
                'Sharp Wave Duration:%5.1f'], durSWR*(1000 / SR))), hold off
            subplot(2, 3, 3), hist(swDiff_chk, 100), ...
                line(repmat(MswDiff+thresSDswD(2)*SDswDiff, [1, 2]), get(gca, 'ylim'), ...
                'linestyle', '--', 'color', 'k')
            line(repmat(swDiff_chk(bound+1), [1, 2]), get(gca, 'ylim'), ...
                'linestyle', '--', 'color', 'r')
            xlabel('A.U.'),
            ylabel('Counts')
            title(sprintf(['Histogram of Sharp Wave Difference Magnitude\n', ...
                '%3.1f SD threshold (black) and detected value (red)'], thresSDswD(2)))
            hold off

            h2 = subplot(2, 3, [4, 5]);
            plot(timeline, ripPow_chk(:, 2)), axis tight, hold on
            plot([timeline(startSWR), timeline(bound+stopSWR)], ...
                [ripPow_chk(startSWR, 2), ripPow_chk(bound+stopSWR, 2)], 'og')
            line([0, 0], get(gca, 'ylim'), ...
                'linestyle', '--', 'color', 'k')
            line(get(gca, 'xlim'), repmat(MripPow, [1, 2]), ...
                'linestyle', '--', 'color', 'k')
            line(get(gca, 'xlim'), repmat(MripPow+thresSDrip(2)*SDripPow, [1, 2]), ...
                'linestyle', '--', 'color', 'r')
            line(get(gca, 'xlim'), repmat(MripPow-thresSDrip(2)*SDripPow, [1, 2]), ...
                'linestyle', '--', 'color', 'r')
            xlabel('Samples'), ylabel('A.U.'), title('Rip Pow'), hold off
            subplot(2, 3, 6), hist(ripPow_chk(:, 2), 100), ...
                line(repmat(MripPow+thresSDrip(2)*SDripPow, [1, 2]), get(gca, 'ylim'), ...
                'linestyle', '--', 'color', 'k')
            line(repmat(ripPow_chk(bound+1, 2), [1, 2]), get(gca, 'ylim'), ...
                'linestyle', '--', 'color', 'r')
            xlabel('A.U.'),
            ylabel('Counts')
            title(sprintf(['Histogram of Ripple Power\n', ...
                '%3.1f SD threshold (black) and detected value (red)'], thresSDrip(2)))
            hold off

            linkaxes([h1, h2], 'x')

            %         keyboard
            pause(.1)

            clf
        end

    end
    SWR_valid.Ts = SWR_valid.Ts(1:valid, :);
    SWR_valid.SwMax = SWR_valid.SwMax(1:valid, :);
    SWR_valid.RipMax = SWR_valid.RipMax(1:valid, :);

    if FIGS
        % Evaluate fixed point Precision/Recall performance
        figure('Color', 'w'),

        % marginal histograms
        [n1, ctr1] = hist(swDiffAll, 1000);
        [n2, ctr2] = hist(ripPowerAll, 1000);

        h1 = subplot(2, 2, 2);
        plot(swDiffAll(idx1), ripPowerAll(idx1), '.r'), hold on
        plot(swDiffAll(idx2), ripPowerAll(idx2), '.b'),
        plot(swDiffAll(SWRindF), ripPowerAll(SWRindF), 'og', 'markerfacecolor', 'none')
        axis([min(swDiffAll(:)), max(swDiffAll(:)), min(ripPowerAll(:)), max(ripPowerAll(:))])
        line(get(h1, 'xlim'), repmat(thresL_Rip, [1, 2]), 'linestyle', '--', 'color', 'k')
        line(repmat(thresL_swD, [1, 2]), get(h1, 'ylim'), 'linestyle', '--', 'color', 'k')
        title('Sharp Wave Difference vs. Ripple Power', 'fontsize', 12)
        xlabel('Sharp Wave Diff', 'fontsize', 12)
        ylabel('Ripple Power', 'fontsize', 12)
        h2 = subplot(2, 2, 4);
        bar(ctr1, -n1, 1),
        axis([min(swDiffAll(:)), max(swDiffAll(:)), -max(n1) * 1.1, 0]),
        axis('off'),
        h3 = subplot(2, 2, 1);
        barh(ctr2, -n2, 1),
        axis([-max(n2) * 1.1, 0, min(ripPowerAll(:)), max(ripPowerAll(:))])
        axis('off'),
        set(h1, 'Position', [0.25, 0.35, 0.55, 0.55], 'fontsize', 12);
        set(h2, 'Position', [.25, .1, .55, .15], 'fontsize', 12);
        set(h3, 'Position', [.05, .35, .15, .55], 'fontsize', 12);

        string = sprintf(['K-means Clustering:\n\n', ...
            '    Center1: X  = %7.2f, Y = %7.2f\n', ...
            '    Center2: X  = %7.2f, Y = %7.2f\n\n', ...
            'Number of Detected SWR:%5d\n'], ...
            Ctrs(1, 1), Ctrs(1, 2), Ctrs(2, 1), Ctrs(2, 2), valid);
        annotation('textbox', [0.81, 0.70, .1, .1], 'String', string, ...
            'edgecolor', 'none', 'fontsize', 12)
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% TIME MINIMUM TROUGH %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Filter pyramidal lfp in ripple band and get the envelope
    lfpfilt = bz_Filter(lfp(:, 1), 'passband', ripBP, 'filter', 'fir1', 'nyquist', SR/2);
    signEnv = 2 * lfpfilt .* lfpfilt;
    envelope = sgolayfilt(signEnv, 4, 101);
    envelope = abs(sqrt(envelope));
    % Compute timestap of the trough nearest to the envelope maximum
    % That timestamp will be the middle timestamp of the ripple output
    TsTrough = SWR_valid.Ts(:, 1) * 0;
    for irip = 1:size(SWR_valid.Ts, 1)
        idxs = SWR_valid.Ts(irip, 2):SWR_valid.Ts(irip, 3);
        [~, imax] = max(envelope(idxs));
        imax = imax + idxs(1) - 1;
        idxsmax = round([imax - 0.01 * SR:imax + 0.01 * SR]);
        [~, imin] = min(lfp(idxsmax));
        imin = imin + idxsmax(1) - 1;
        TsTrough(irip) = imin;
        %     plot((idxs-imin)/SR, lfp(idxs,1), 'k', 'linewidth', 2)
        %     hold on
        %     plot((idxs-imin)/SR, lfpfilt(idxs))
        %     plot((idxs-imin)/SR, 1500*signEnv(idxs)/max(signEnv(idxs)))
        %     plot(0/SR, lfp(imin), '*')
        %     hold off
        %     xlim([-0.05 0.05])
        %     pause
    end

    % Keep ripples whose computed troughs are inside the start and end timestamps
    idxsIn = find(TsTrough > SWR_valid.Ts(:, 2) & TsTrough < SWR_valid.Ts(:, 3));

    % Compute instantaneous phase and amplitude
    h = hilbert(lfpfilt);
    phase = angle(h);
    amplitude = abs(h);
    unwrapped = unwrap(phase);
    % Compute instantaneous frequency
    frequency = bz_Diff(medfilt1(unwrapped, 12), [0:length(lfpfilt) - 1]'/SR, 'smooth', 0);
    frequency = frequency / (2 * pi);

    %%%%%%%%%%%%%%%%%%%%%%%
    %%% GENERATE OUTPUT %%%
    %%%%%%%%%%%%%%%%%%%%%%%
    % 1) swr detections
    % - Start and end timestamps of each SWR
    SWR.timestamps = SWR_valid.Ts(idxsIn, [2, 3]) / SR;
    % Duration
    SWR.duration = SWR.timestamps(:, 2) - SWR.timestamps(:, 1);
    % amplitude at trough
    SWR.amplitude = amplitude(TsTrough(idxsIn));
    % Frequency at trough
    SWR.frequency = frequency(TsTrough(idxsIn));
    % - Peak power of ripple timestamp of each SWR
    SWR.peaks = TsTrough(idxsIn) / SR; %SWR_valid.Ts(:,1)/SR;
    % - Peak power values of each SWR
    SWR.peakNormedPower = [];
    % - Name of this function
    SWR.detectorName = 'DetectSWR';
    % - Standard deviation used to threshold ripple power
    SWR.stdev = [thresSDrip(2), thresSDswD(2)];
    SWR.SwMax = SWR_valid.SwMax(idxsIn);
    SWR.RipMax = SWR_valid.RipMax(idxsIn);
    SWR_valid.Ts = SWR_valid.Ts(idxsIn, :);
    SWR_valid.SwMax = SWR_valid.SwMax(idxsIn, :);
    SWR_valid.RipMax = SWR_valid.RipMax(idxsIn, :);

    % Make sure all ripple epochs are unique.
    % check start ts
    [~, unique_idx, ~] = unique(SWR.timestamps(:, 1));
    SWR.timestamps = SWR.timestamps(unique_idx, :);
    SWR.duration = SWR.duration(unique_idx, :);
    SWR.amplitude = SWR.amplitude(unique_idx, :);
    SWR.frequency = SWR.frequency(unique_idx, :);
    SWR.peaks = SWR.peaks(unique_idx, :);
    SWR.SwMax = SWR.SwMax(unique_idx, :);
    SWR.RipMax = SWR.RipMax(unique_idx, :);
    SWR_valid.Ts = SWR_valid.Ts(unique_idx, :);
    SWR_valid.SwMax = SWR_valid.SwMax(unique_idx, :);
    SWR_valid.RipMax = SWR_valid.RipMax(unique_idx, :);
    % check peak ts
    [~, unique_idx, ~] = unique(SWR.peaks);
    SWR.timestamps = SWR.timestamps(unique_idx, :);
    SWR.duration = SWR.duration(unique_idx, :);
    SWR.amplitude = SWR.amplitude(unique_idx, :);
    SWR.frequency = SWR.frequency(unique_idx, :);
    SWR.peaks = SWR.peaks(unique_idx, :);
    SWR.SwMax = SWR.SwMax(unique_idx, :);
    SWR.RipMax = SWR.RipMax(unique_idx, :);
    SWR_valid.Ts = SWR_valid.Ts(unique_idx, :);
    SWR_valid.SwMax = SWR_valid.SwMax(unique_idx, :);
    SWR_valid.RipMax = SWR_valid.RipMax(unique_idx, :);

    try
    %Write event file: Generate Output and Write out Event file
    if EVENTFILE
        rippleFiles = dir('*.R*.evt');
        if isempty(rippleFiles)
            fileN = 1;
        else
            %set file index to next available value\
            pat = '.R[0-9].';
            fileN = 0;
            for ii = 1:length(rippleFiles)
                token = regexp(rippleFiles(ii).name, pat);
                val = str2double(rippleFiles(ii).name(token+2:token+4));
                fileN = max([fileN, val]);
            end
            fileN = fileN + 1;
        end
        fid = fopen(sprintf('%s%s%s.R%02d.evt', pathname, filesep, filename, fileN), 'w');

        % convert detections to milliseconds
        SWR_final = SWR_valid.Ts * (1000 / SR);
        fprintf(1, 'Writing event file ...\n');
        for ii = 1:size(SWR_final, 1)
            fprintf(fid, '%9.1f\tstart\n', SWR_final(ii, 2));
            fprintf(fid, '%9.1f\tpeak\n', SWR_final(ii, 1));
            fprintf(fid, '%9.1f\tstop\n', SWR_final(ii, 3));
        end
        fclose(fid);
    end
    end

    % 3) analysis params
    params.filename = filename;
    params.Channels = Channels;
    params.SR = SR;
    params.nChan = nChan;
    params.swBP = swBP;
    params.ripBP = ripBP;
    params.per_thresswD = per_thresswD;
    params.per_thresRip = per_thresRip;
    params.WinSize = WinSize;
    params.Ns_chk = Ns_chk;
    params.thresSDswD = thresSDswD;
    params.thresSDrip = thresSDrip;
    params.minIsi = minIsi;
    params.minDurSW = minDurSW;
    params.maxDurSW = maxDurSW;
    params.minDurRP = minDurRP;
    params.EVENTFILE = EVENTFILE;
    params.TRAINING = TRAINING;
    params.DEBUG = DEBUG;
    SWR.detectorinfo.detectionparms = params;
    SWR.detectorinfo.detectorname = 'bz_DetectSWR';
    SWR.detectorinfo.detectiondate = datestr(datetime('now'), 'dd-mmm-yyyy');
    SWR.detectorinfo.detectionintervals = Epochs;

    try
    % write out log file
    logAnalysisFile(mfilename('fullpath'), pathname);
    end
    
    ripples = SWR;

    try
    if saveMat
        save([Filebase, '.ripples.events.mat'], 'ripples')
    end
    end
end


% ------------------------------- Helper functions -------------------------------

function gwin = makegausslpfir(Fc, Fs, s)

% makegausslpfir        Gaussian LP filter
%
% win = makegausslpfir( Fc, Fs, s )
%
% Fc    corner frequency [Hz]
% Fs    sampling frequency [Hz]
% s     support [SD], minimum 3, default 4
%
% win   fir
%
% see also  firfilter
% 02-oct-13 ES

nargs = nargin;
if nargs < 1 || isempty(Fc), Fc = 100; end
if nargs < 2 || isempty(Fs), Fs = 1000; end
if nargs < 3 || isempty(s), s = 4; end
s = max(s, 3);

sd = Fs / (2 * pi * Fc);
x = -ceil(s*sd):ceil(s*sd);
gwin = 1 / (2 * pi * sd) * exp(-(x.^2 / 2 / sd.^2));
gwin = gwin / sum(gwin);


% EOF

% to determine the proper support:
% cwin = cumsum( gwin );
% length( cwin ), [ find( cwin >= 0.001, 1, 'first' ) - 1, find( cwin <= 0.999, 1, 'last' ) + 1 ]
% 3 SD are usually enough


% FIRFILT       FIR filtering with zero phase distortion.
%
%               matrix columns are filtered;
%               output is a column vector / matrix.
% 19-Jul-02 ES
% 27-jan-03 zero phase lag

function Y = firfilt(x, W)

if all(size(W) >= 2), error('window must be a vector'), end
if numel(x) == max(size(x)), x = x(:); end

C = length(W);
if C > size(x, 1)
    Y = NaN * ones(size(x));
    return
end
D = ceil(C/2) - 1;
Y = filter(W, 1, [flipud(x(1:C, :)); x; flipud(x(end-C+1:end, :))]);
clear x
Y = Y(1+C+D:end-C+D, :);



function logAnalysisFile(AnalysisFileName, writePath)
% The purpose of this utility is to provide a record of analyses run by
% writing out to a log file the state of the code that was used to analyze
% the data. It writes the log file into the current directory if only a
% filename is provided. An absolute path is required to write out to an
% arbitrary location.
%
% This utility is designed to be placed within an .m file used for
% analysis.
%
% dependencies: mfilename
%
% author: John D. Long II, PhD   contact: jlong29@gmail.com
%
%%%%%%%%%%%%%%
%%% INPUTS %%%
%%%%%%%%%%%%%%
% AnalysisFileName: the output of mfilename('fullname') called from within
%   the .m file
% writePath (optional): a user specified path for writing the log file.
%
%%%%%%%%%%%%%%%
%%% OUTPUTS %%%
%%%%%%%%%%%%%%%
% A .log file entitled mfilename-date.log written to the location
%
% Example usage:
%   logAnalysisFile(mfilename('fullpath') or
%   logAnalysisFile(mfilename('fullpath','~/mnt/Analysis/')

% input check
if nargin < 2 || isempty(writePath)
    writePath = [];
end

% assumes output from mfilename
info         = dir([AnalysisFileName '.m']);
[~, filename] = fileparts([AnalysisFileName, '.m']);

% Check for log file and create or update as required
if ~isempty(writePath)
    % check if writePath is a valid directory
    if exist(writePath, 'dir')
        % check if writePath ends in filesep
        if ~strcmp(writePath(end), filesep)
            writePath = [writePath, filesep];
        end
        % check for pre-existing log file in writePath
        loginfo = dir(sprintf('%s%s.log', writePath, filename));
        if isempty(loginfo)
            fprintf(1, 'First Use Of Analysis File: Writing Log.\n');
            fid1 = fopen(sprintf('%s%s.log', writePath, filename), 'w+');
        else
            % let's check if we need to update this file
            fid1 = fopen(sprintf('%s%s.log', writePath, filename), 'r');
            tline = fgetl(fid1);
            if strcmp(tline, info.date)
                fprintf(1, 'Analysis File Up To Date.\n');
                return
            else
                fprintf(1, 'Analysis File Changed: Writing Log.\n');
                fclose(fid1);
                fid1 = fopen(sprintf('%s%s.log', writePath, filename), 'w+');
            end
        end
    else

        fprintf(1, 'InputWarning: Log File Written out to launch directory.\n');

        % check for pre-existing log file current directory
        loginfo = dir(sprintf('%s.log', filename));
        if isempty(loginfo)
            fprintf(1, 'First Use Of Analysis File: Writing Log.\n');
            fid1 = fopen(sprintf('%s.log', filename), 'w+');
        else
            % let's check if we need to update this file
            fid1 = fopen(sprintf('%s.log', filename), 'r');
            tline = fgetl(fid1);
            if strcmp(tline, info.date)
                fprintf(1, 'Analysis File Up To Date.\n');
                return
            else
                fprintf(1, 'Analysis File Changed: Writing Log.\n');
                fclose(fid1);
                fid1 = fopen(sprintf('%s.log', filename), 'w+');
            end
        end
    end
else
    % check for pre-existing log file current directory
    loginfo = dir(sprintf('%s.log', filename));
    if isempty(loginfo)
        fprintf(1, 'First Use Of Analysis File: Writing Log.\n');
        fid1 = fopen(sprintf('%s.log', filename), 'w+');
    else
        % let's check if we need to update this file
        fid1 = fopen(sprintf('%s.log', filename), 'r');
        tline = fgetl(fid1);
        if strcmp(tline, info.date)
            fprintf(1, 'Analysis File Up To Date.\n');
        else
            fprintf(1, 'Analysis File Changed: Writing Log.\n');
            fclose(fid1);
            fid1 = fopen(sprintf('%s.log', filename), 'w+');
        end
    end
end
% Open the analysis file to be read
fid2 = fopen(sprintf('%s', [AnalysisFileName, '.m']), 'r');
fprintf(fid1, '%s\n', info.date);

while ~feof(fid2)
    try
        tline = fgetl(fid2);
        fprintf(fid1, '%s\n', tline);
    catch
        keyboard
    end
end
% Close all files
fclose(fid1);
fclose(fid2);


function Y = FiltFiltM(b, a, X, Dim)
% FiltFiltM - zero phase-shift digital filter
% The signal is filtered in forward and reverse direction for a zero phase
% distortion. This is equivalent to FILTFILT of Matlab's Signal Processing
% Toolbox, but the direct processing of arrays let this function run 10% to 90%
% faster. With the fast Mex function FilterX, FiltFiltM is 80% to 95% faster
% than FILTFILT and uses much less temporary memory. 
%
% Y = FiltFiltM(b, a, X, Dim)
% INPUT:
%   b, a: Filter coefficients as vectors. See FILTER.
%   X:    Signal as DOUBLE or SINGLE array.
%   Dim:  Dimension to operate on. If Dim is the empty matrix [], the filter is
%         applied to the 1st non-singelton dimension. Filtering along the 1st
%         dimension is faster than for the other dimensions.
%         Optional, default: 1.
% OUTPUT:
%   Y:    Array of same size and type as X. The values are filtered in forward
%         and backward direction for a zero phase shift.
%         The result equals the output of FILTFILT.
%
% NOTES:
% - In opposite to FILTFILT the initial and final phases are treated separately
%   instead of enlarging the signal array. This reduces the memory consumption.
% - The results for signals of type SINGLE are more accurate, because FilterX
%   uses DOUBLE values for computing the filter.
%
% EXAMPLES:
%   b = [0.0675, 0.1349, 0.0675];  % Low pass Butterworth filter
%   a = [1.0,   -1.1430, 0.4128];
%   X = rand(16, 8);
%   Y = FiltFiltM(b, a, X);        % Along the 1st dimension
%   figure, subplot(1,2,1); plot(Y);
%   Y = FiltFiltM(b, a, X, 2);     % Along the 2nd dimension
%   subplot(1,2,2), plot(Y.');
%
% INSTALLATION:
%   This function runs without installation. But for more speed the C-mex
%   FilterX.c can be compiled: Call FilterM() without inputs. See there for
%   further instructions.
%
% See also FILTER, FILTER2, FILTFILT, CONV, FilterM.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, 2009-2011 matlab.THISYEAR(a)nMINUSsimon.de

% References:
% [1] Sanjit K. Mitra, Digital Signal Processing, 2nd ed., McGraw-Hill, 2001
% [2] Fredrik Gustafsson, Determining the initial states in forward-backward
%     filtering, IEEE Transactions on Signal Processing, pp. 988--992,
%     April 1996, Volume 44, Issue 4

% $JRev: R-t V:019 Sum:flRzmua2/Vps Date:20-Jul-2011 10:29:21 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_FiltFiltM $
% $File: Tools\GLMath\FiltFiltM.m $
% History:
% 009: 08-May-2011 16:45, Get initial and final parts separately: 10% faster.
%      Joining the reflected initial and final parts to the signal wastes
%      memory. Now the reflected edge parts are filtered separately just to get
%      the initial conditions for filtering the original signal.
% 011: 12-Jul-2011 16:40, 4th input [Dim]. Rewritten from scratch.
%      Now this function does not contain any parts of the original FILTFILT
%      anymore. The initial conditions are determined following Gustafsson's
%      article.

% Initialize: ==================================================================
% Global Interface: ------------------------------------------------------------
% Use the much faster C-Mex replacement for FILTER:
persistent hasFilterX
if isempty(hasFilterX)
   hasFilterX = ~isempty(which('FilterX'));
   if ~hasFilterX && ~isempty(which('FilterM'));
      fprintf(['::: %s: Matlab''s FILTER is used.\n', ...
         '    Run FilterM without inputs to compile FilterX.c\n'], ...
         mfilename);
   end
end

% Initial values: --------------------------------------------------------------
% Program Interface: -----------------------------------------------------------
nArg = nargin;
if nArg == 3
   Dim = [];    % Default is specified later: 1st or 1st non-singelton dimension
elseif nArg == 4
   if numel(Dim) > 1 || ~isequal(Dim, round(abs(Dim)))
      error(['JSimon:', mfilename, ':BadDim'], ...
         ['*** ', mfilename, ': Dim must be a valid dimension.']);
   end
else
   error(['JSimon:', mfilename, ':BadNInput'], ...
      ['*** ', mfilename, ': 3 or 4 inputs required.']);
end

if isempty(b) || isempty(a) || isempty(X)
   Y = X([]);
   return;
end

% Get filter parameter vectors as row vectors and padded with zeros on demand:
nb    = numel(b);
na    = numel(a);
b     = b(:);
a     = a(:);
Order = max(nb, na);
if nb < Order
   b(Order, 1) = 0;
end
if na < Order
   a(Order, 1) = 0;
end
nEdge = 3 * (Order - 1);   % Number of edge elements

% Move dimension to operate on to dim 1: ---------------------------------------
% Specifying the dimension in the FILTER command is not faster. It seems like
% the permutation is performed inside FILTER also instead of operating on the
% original data.
sizeX  = size(X);
ndimsX = ndims(X);
if ndimsX == 2             % Vector or matrix
   if isempty(Dim)         % Default: First non-singelton dimension
      doPermute = false;
      doReshape = (sizeX(1) == 1);
      if doReshape
         X = X(:);
      end
   elseif Dim < 3
      doReshape = false;
      doPermute = (Dim == 2);
      if doPermute
         perm = [2, 1];
         X    = transpose(X);
      end
   else
      % Usually Matlab functions allow appending arbitrary trailing singelton
      % dimensions, but then the (xLen<=nEdge) check fails anyway:
      error(['JSimon:', mfilename, ':BadDim'], ...
         ['*** ', mfilename, ': Dim out of range.']);
   end
   
else                % Signal has more than 2 dimensions: -----------------------
   doReshape = true;
   if isempty(Dim)         % First non-singleton dimension:
      doPermute = false;
      nonSingle = find(sizeX ~= 1);
      if ~isempty(nonSingle) && nonSingle(1) > 1  % Short-circuit needed
         X = reshape(X, sizeX(nonSingle(1)), []);
      end
   elseif Dim <= ndimsX    % Dim specified:
      if Dim == 1
         doPermute = false;
      else
         doPermute = true;
         perm      = [Dim, 1:Dim - 1, Dim + 1:ndims(X)];
         X         = reshape(permute(X, perm), sizeX(Dim), []);
         sizeX     = sizeX(perm);
      end
   else
      error(['JSimon:', mfilename, ':BadDim'], ...
         ['*** ', mfilename, ': Dim out of range.']);
   end
end

xLen = size(X, 1);          % Not sizeX - would fail for row vectors!
if xLen <= nEdge            % Input data too short
   error(['JSimon:', mfilename, ':BadDataLen'], ...
      ['*** ', mfilename, ...
      ': Signal must be longer than 3 times the filter order.']);
end

% User Interface: --------------------------------------------------------------
% Do the work: =================================================================
% Create initial conditions to treat offsets at beginning and end:
K       = eye(Order - 1);
K(:, 1) = a(2:Order);
K(1)    = K(1) + 1;
K(Order:Order:numel(K)) = -1.0;
IC      = K \ (b(2:Order) - a(2:Order) * b(1));

% Use a reflection to extrapolate signal at beginning and end to reduce edge
% effects (BSXFUN would be some micro-seconds faster, but it is not available in
% Matlab 6.5):
x1_2 = 2 * X(1,    :);
xf_2 = 2 * X(xLen, :);
Xi   = x1_2(ones(1, nEdge), :) - X((nEdge + 1):-1:2, :);
Xf   = xf_2(ones(1, nEdge), :) - X((xLen - 1):-1:(xLen - nEdge), :);

if hasFilterX  % Use the faster C-mex filter function: -------------------------
   % Filter initial reflected signal:
   [dum, Zi] = FilterX(b, a, Xi, IC * Xi(1, :));   %#ok<ASGLU>
   
   % Use the final conditions of the initial part for the actual signal:
   [Ys, Zs]  = FilterX(b, a, X,  Zi);              % "s"teady state
   Yf        = FilterX(b, a, Xf, Zs);              % "f"inal conditions
   
   % Filter signal again in reverse order:
   [dum, Zf] = FilterX(b, a, Yf, IC * Yf(nEdge, :), 'reverse');  %#ok<ASGLU>
   Y         = FilterX(b, a, Ys, Zf, 'reverse');
   
else  % Use the slower built-in FILTER function: -------------------------------
   [dum, Zi] = filter(b, a, Xi, IC * Xi(1, :));     %#ok<ASGLU>
   [Ys, Zs]  = filter(b, a, X,  Zi);                % "s"teady state
   Yf        = filter(b, a, Xf, Zs);                % "f"inal conditions
   Yf        = Yf(nEdge:-1:1, :);
   [dum, Zf] = filter(b, a, Yf, IC * Yf(1, :));     %#ok<ASGLU>
   Y         = filter(b, a, Ys(xLen:-1:1, :), Zf);  % Filter reverted signal
   Y         = Y(xLen:-1:1, :);                     % Re-revert
end

% Reconstruct original array dimension: ----------------------------------------
if doReshape
   Y = reshape(Y, sizeX);
end
if doPermute
   Y = ipermute(Y, perm);
end


% return;

function filtered = bz_Filter(samples,varargin)

%Filter - bandpass filter samples
%
%  USAGE
%
%    filtered = bz_Filter(samples,<options>)
%
%INPUT
%    samples        samples given as a [t x channels] timeseries vector
%                   or as buzcode structure with fields samples.data,
%                   samples.timestamps, samples.samplingRate
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'passband'    pass frequency range. [0 X] for low-pass, [X inf] for highpass
%     'stopband'    stop frequency range
%     'order'       filter order (number of cycles, default = 4)
%     'ripple'      filter ripple (default = 20)
%     'filter'      choose filter type between 'cheby2' (default) and 'fir1'
%     'nyquist'     nyquist frequency (default = 625), calculated
%                   automatically from samples.samplingRate for BUZCODE
%     'FMAlegacy'   true/false, uses FMA legacy input style
%                   (samples given as a list of (t,v1,v2,v3...) tuples)
%     'fast'        true/false, uses FiltFiltM if true (doesn't work w/ new
%                   version of matlab)
%     'channels'    if input is a buzcode lfp structure with field
%                   samples.channels, will only filter the selected
%                   channels
%     'intervals'   only filter in given intervals (only works for buzcode input)
%    =========================================================================
%
%OUTPUT
%   filtered        -if the input is a timeseries vector, output is as well
%                   -if the input is a buzcode structure, output is a
%                    buzcode structure with the following fields:
%       .data       the filtered data
%       .phase      phase, calculated by the hilbert method
%       .amp        amplitude, calculated by the hilbert method
%       .timestamps     (s)
%       .samplingRate   (Hz)
%       .filterparms    a structure that holds the parameters used for
%                       filtering, for future reference
%
%   Dependencies: InIntervals, FiltFiltM
%
% Copyright (C) 2004-2011 by Michaël Zugaro
% updated 2017 DLevenstein for buzcode
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%% Parms and Defaults

% Default values
passband = [];
stopband = [];
order = 4;
ripple = 20;
nyquist = 625; %written over later from samples.samplingRate if available
type = 'cheby2';
FMAlegacy = false;
BUZCODE = false;
fast = false;
channels = [];
intervals = [-Inf Inf];

% Check number of parameters
if nargin < 1 || mod(length(varargin),2) ~= 0
    error('Incorrect number of parameters (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
end

if isempty(samples)
    filtered = samples;
    return;
end

% Parse parameter list
for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i) ' is not a property (type ''help <a href="matlab:help Filter">Filter</a>'' for details).']);
    end
    switch(lower(varargin{i}))
        case 'passband'
            if ~isempty(stopband)
                error('Cannot specify both a passband and stopband (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
            end
            passband = varargin{i+1};
            if ~isdvector(passband,'#2','>=0')
                error('Incorrect value for ''passband'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
            end
        case 'stopband'
            if ~isempty(passband)
                error('Cannot specify both a passband and stopband (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
            end
            stopband = varargin{i+1};
            if ~isdvector(stopband,'#2','>=0')
                error('Incorrect value for ''stopband'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
            end
        case 'filter'
            type = lower(varargin{i+1});
            if ~isastring(type,'cheby2','fir1','butter')
                error(['Unknown filter type ''' type ''' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).']);
            end
            
        case 'order'
            order = lower(varargin{i+1});
            if ~isiscalar(order,'>0')
                error('Incorrect value for ''order'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
            end
            
        case 'ripple'
            ripple = lower(varargin{i+1});
            if ~isiscalar(ripple,'>0')
                error('Incorrect value for ''ripple'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
            end
            
        case 'nyquist'
            nyquist = varargin{i+1};
            if ~isiscalar(nyquist,'>0')
                error('Incorrect value for property ''nyquist'' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
            end
        case 'fast'
            fast = varargin{i+1};
            if ~islogical(fast)
                error('Incorrect value for property ''FMALegacy''');
            end
        case 'fmalegacy'
            FMAlegacy = varargin{i+1};
            if ~islogical(FMAlegacy)
                error('Incorrect value for property ''FMALegacy''');
            end
            
        case 'channels'
            channels = varargin{i+1};
            
        case 'intervals'
            intervals = varargin{i+1};
            
        otherwise
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help Filter">Filter</a>'' for details).']);
            
    end
end
%%
%Check the input structure for buzcode stuff
if isstruct(samples)
    BUZCODE = true;
    if ~isfield(samples,'data')
        error('Your input structure has no .data field :(')
    end
    if ~isfield(samples,'samplingRate')
        warning('Your input strucutre has no .samplingRate field, using provided/default nyquist frequency')
        samples.samplingRate = 2.*nyquist;
    end
    if ~isfield(samples,'timestamps')
        warning('Your input structure has no .timestamps field, one will be provided for you')
        samples.timestamps = (0:length(samples.data)-1)./samples.samplingRate;
    end
    
    if length(samples.samplingRate) > 1
        nyquist = 0.5.*samples.samplingRate(1);
    else
        nyquist = 0.5.*samples.samplingRate;
    end
end

if isempty(passband) && isempty(stopband)
    error('Missing passband or stopband (type ''help <a href="matlab:help Filter">Filter</a>'' for details).');
end


switch(type)
    case 'cheby2'
        if ~isempty(passband)
            if passband(1) == 0
                [b, a] = cheby2(order,ripple,passband(2)/nyquist,'low');
            elseif passband(2) == inf
                [b, a] = cheby2(order,ripple,passband(1)/nyquist,'high');
            else
                [b, a] = cheby2(order,ripple,passband/nyquist);
            end
        else
            [b, a] = cheby2(order,ripple,stopband/nyquist,'stop');
        end
        warning('Cheby2 is often numerically unstable - if you get NaNs, use the ''fir1'' argument to bz_Filter');
    case 'fir1'
        %Order input to fir1 needs to be in samples, not cycles
        if ~isempty(passband)
            if passband(1) == 0
                filt_order = round(order*2*nyquist./passband(2));
                [b, a] = fir1(filt_order,passband(2)/nyquist,'low');
            elseif passband(2) == inf
                filt_order = round(order*2*nyquist./passband(1));
                [b, a] = fir1(filt_order,passband(1)/nyquist,'high');
            else
                filt_order = round(order*2*nyquist./passband(1));
                [b, a] = fir1(filt_order,passband/nyquist);
            end
        else
            filt_order = round(order*2*nyquist./stopband(1));
            [b, a] = fir1(filt_order,stopband/nyquist,'stop');
        end
    case 'butter'
        if ~isempty(passband)
            [b, a] = butter(order,[passband(1)/nyquist passband(2)/nyquist],'bandpass');
        else
            [b, a] = butter(order,stopband(1)/nyquist,'stop');
        end
end


if FMAlegacy %FMA has (samples given as a list of (t,v1,v2,v3...) tuples)
    filtered(:,1) = samples(:,1);
    for i = 2:size(samples,2)
        if ~fast
            filtered(:,i) = filtfilt(b,a,double(samples(:,i)));
        else
            filtered(:,i) = FiltFiltM(b,a,double(samples(:,i)));
        end
    end
    
elseif BUZCODE %BUZCODE has samples as a data structure
    
    %Restrict to intervals, with overhang to remove edge effects at transitions
    %(then remove later)
    overhang = (order)./passband(1);
    overint = bsxfun(@(X,Y) X+Y,intervals,overhang.*[-1 1]);
    keepIDX = InIntervals(samples.timestamps,overint);
    samples.data = samples.data(keepIDX,:);
    filtered.timestamps = samples.timestamps(keepIDX);
    
    
    if ~isempty(channels)
        [~,usechannels] = ismember(channels,samples.channels);
        samples.data = samples.data(:,usechannels);
        filtered.channels = samples.channels(usechannels);
    end
    for i = 1:size(samples.data,2)
        if ~fast
            filtered.data(:,i) = filtfilt(b,a,double(samples.data(:,i)));
        else
            filtered.data(:,i) = FiltFiltM(b,a,double(samples.data(:,i)));
        end
        hilb = hilbert(filtered.data(:,i));
        filtered.hilb(:,i) = hilb;
        filtered.amp(:,i) = abs(hilb);
        filtered.phase(:,i) = angle(hilb);
    end
    
    %Remove the overhang from intervals
    keepIDX = InIntervals(filtered.timestamps,intervals);
    filtered.data = filtered.data(keepIDX,:);
    filtered.hilb = filtered.hilb(keepIDX,:);
    filtered.amp = filtered.amp(keepIDX,:);
    filtered.phase = filtered.phase(keepIDX,:);
    filtered.timestamps = filtered.timestamps(keepIDX);
    
    filtered.filterparms.passband = passband;
    filtered.filterparms.stopband = stopband;
    filtered.filterparms.order = order;
    filtered.samplingRate = samples.samplingRate;
else %or if you just want filter a basic timeseries
    for i = 1:size(samples,2)
        if ~fast
            filtered(:,i) = filtfilt(b,a,double(samples(:,i)));
        else
            filtered(:,i) = FiltFiltM(b,a,double(samples(:,i)));
        end
        
    end
end


function basename = basenameFromBasepath(basepath)
% basenameFromBasepath-gets basename from a basepath input to a function.  Uses name of the last
% directory on the path to give basename.  Takes care of the problem of the
% fact that if a fileseparator is on the end of the basepath, you have to
% take it off first.
%
% INPUTS
% basepath  is assumed to be a directory containing all data for a given 
%           recording session.  It is assumed organization is 
%           /whateverdirectory/basename/basename.xml and basepath is 
%           equivalent to whateverdirectory/basename/.  
%
% OUTPUTS
% basename  as specified above
%
% Brendon Watson, 2017


if strcmp(basepath(end),filesep);
    basepath = basepath(1:end-1);
end
% if there are dots in the filename take that last folder
if contains(basepath,'.')
    parts = strsplit(basepath,filesep);
    basename = parts{end};
    return
end
[~,basename] = fileparts(basepath);



function z = xmltools( arg, out_file, varargin)
% XMLTOOLS - tools for managing xml data sets
%      - if arg is a string : arg is an XML file to convert into MATLAB struct
%      - if arg is a variable : it is a MATLAB struct to write into XML, to stdout if out_file is not given
% use :
%  z = xmltools('filename.xml'); read an xml file and store it into z
%  xmltools(z,'filename.xml'); write z into the file
%  xmltools(z,'get','tag-name'); returns only subset of z child which name is tag-name
%  xmltools(z,'get-attrib', 'attrib-name')
%  xmltools(z,'get','tag-name', 'attribs'|'value'|'child');
%
% project 'File parsing'
% title    'XML parsing'
% author  'Charles-Albert Lehalle'
% mailto  'charles.lehalle@miriadtech.com'
% version '2.5'
% date    'mar2003--sept2003'

version = '2.5';

%%** XML TOOLS FOR MATLAB
% This is an OCAMAWEB (http://ocamaweb.sourceforge.net) generated documentation.\\
% This function manage the exchange of data between XML files and a MATLAB structure.
% The MATLAB structure is this : a {\bf node} is a struct with fields {\it tag}, {\it value}, {\it attribs}, {\it child}.
% Where {\it tag} is the tag name, {\it value} its contents, {\it attribs} an array of structure with fields {\it name} and {\it value},
% and {\it child} an array of such nodes.\\
% All those fields are always present and can be empty except for the first node (root) that has only children.\\
%

global Verbose;
Verbose=0;
%%* READ AN XML FILE
if isstr(arg)
  
  %< R�cup�ration du fichier dans un string
  fid = fopen(arg, 'r');
  if fid<0
      error('Pb opening the xml file');
  end
  F = fread(fid);
  s = char(F');
  fclose(fid);
  %>
  
  %< Parsing
  z = parse_xml(s);
  %>
  return
end

if ~isstr(arg)

  %<* SELECT A SUBSET OF z CHILD
  if length(varargin) >= 1 %& ~FileExists(arg)
    
    cmode = upper( out_file);
    z     = arg;
    
    switch upper( cmode)
      
     case 'GET'
      %< Get the subset
      % warning: I will have to change the value of next at some places
      next = 'child';
    
      if ~isfield(z, next)
        error('XMLTOOLS:GET', 'For child selection, structured first argument is needed');
      end
      tag_name = (varargin{1});
      
      z = get_childs(z, next, tag_name);
      %>
      
      %< get values
      % xmltools( z, 'get', 'tag-name', 'attribs'|'value'|'child')
      if length(varargin) == 2
        switch upper( varargin{2})
         case 'VALUE'
          z = [z.value];
         case 'ATTRIBS'
          z = [z.attribs];
         case 'CHILD'
          z = [z.child];
         otherwise
          error('XMLTOOLS:GET', 'get mode <%s> is not defined, use one of <attribs>|<value>|<child>', upper(varargin{2}));
        end
      end
      %>
      
     case 'GET-ATTRIB'
      %< get attrib
      % xmltools(z, 'get-attrib', 'attrib-name')
      s = z.attribs;
      for i=1:length(s)
        if strcmp((s(i).name), ( varargin{1}) )
          z = s(i).value;
          return
        end
      end
      error('XMLTOOLS:GET-ATTRIB', 'no attribute found'); %, varargin{1});
      %>
      
    end
           
    return
  end
  %>*
  
  %<* WRITE AN XML STRUCTURE
  
  %< Selection de la cible
  if nargin < 2
    fid = 1;
  else
    fid = fopen(out_file, 'w');end
  end
  %>
  
  %< Ecriture proprement dite
  write_xml(fid, arg, -1);
  %>
  
  %< Fermeture
  if nargin > 1
    fclose(fid);
  end
  %>
  
  %>*


%%** Fonctions internes

%<* parser un string xml
function [z, str] = parse_xml( str, current_tag, current_value, attribs, idx)
global Verbose;
next = 'child';

if nargin < 2
  current_tag   = '';
  current_value = '';
  attribs       = '';
  idx           = 0;
end
z = [];

eot = 0;

while ~eot & ~isempty(udeblank(deblank(str)))
  
  f_end = strfind(str, '</');
  f_beg = strfind(str, '<');
  
  %< Si je n'ai plus de tag dans mon document
  if isempty(f_end) & isempty(f_beg)
    
    if ~strcmp(lower(current_tag), '?xml') & ~isempty(current_tag)
      error('xmltools:parse_xml', 'malformed xml string (current [%s])', current_tag);
    else
      fprintf('end parsing at level %d\n',idx);
      eot = 1;
      return
    end
  end
  %>
  
  if isempty(f_end)
    f_end = length(str)
  else
    f_end = f_end(1);
  end
  if isempty(f_beg)
    f_beg = length(str)
  else
    f_beg = f_beg(1);
  end
  
  if f_end <= f_beg
    %< je rencontre une fermeture
    new_tag = str((f_end+2):end);
    str_t   = str(1:f_end-1);
    f_end = strfind(new_tag,'>');
    if isempty(f_end)
      error('xmltools:parse_xml', 'malformed xml string : never ending tag [%s] encountered', current_tag);
    end
    f_end = f_end(1);
    str     = new_tag(f_end+1:end); % reste
    new_tag = new_tag(1:f_end-1);
    if ~strcmp((new_tag), (current_tag))
      error('xmltools:parse_xml', 'malformed xml string : [%s] not properly closed (closing [%s] encountered)', current_tag, new_tag);
    end
    if Verbose
       fprintf('%sclose [%s]\n', repmat(' ', 2*(idx-1),1), current_tag);
    end
    z.tag     = (current_tag);
    z.attribs = parse_attribs(attribs);
    z.value   = udeblank(deblank(sprintf('%s %s',current_value, str_t)));
    eot       = 1;
    %>
  else
    %< je rencontre une ouverture
    % je vais appeler le m�me code sur ce qu'il y a apr�s moi
    current_value = sprintf('%s %s', current_value, str(1:f_beg-1));
    new_tag   = str(f_beg+1:end);
    f_end = strfind(new_tag,'>');
    if isempty(f_end)
      error('xmltools:parse_xml', 'malformed xml string : never ending tag encountered');
    end
    f_end   = f_end(1);
    str_t   = new_tag(f_end+1:end);
    new_tag = new_tag(1:f_end-1);
    if (new_tag(end) == '/')|(new_tag(end) == '?')
      %< Self closing tag
      % Je met (temporairement!) eot � 1, cela me permet de passer quelques lignes
      % de code tranquilement
      eot = 1;
      %>
    end
    %< Attributs
    f_beg   = strfind(new_tag, ' ');
    if isempty(f_beg)
      new_attribs = '';
      if eot
	new_tag = new_tag(1:end-1);
      end
    else
      new_attribs = new_tag(f_beg+1:end);
      if eot
	new_attribs = new_attribs(1:end-1);
      end
      new_tag     = new_tag(1:f_beg-1);
    end
    %>
    if Verbose
        fprintf('%sopen  [%s]\n', repmat(' ', 2*idx,1), new_tag);
    end
    if eot
      %< If self-colsing tag
      if Verbose
          fprintf('%sclose [%s]\n', repmat(' ', 2*idx,1), new_tag);
      end
      new_attribs = parse_attribs( new_attribs);
      if isfield(z, next)
        nxt = getfield(z, next);
        nxt(end+1) = struct( 'tag', new_tag, 'attribs', new_attribs, 'value', '', next, []);
        z   = setfield(z, next, nxt);
	%z.(next)(end+1) = struct( 'tag', new_tag, 'attribs', new_attribs, 'value', '', next, []);
      else
        z = setfield(z, next, struct( 'tag', new_tag, 'attribs', new_attribs, 'value', '', next, []) );
	%z.(next) = struct( 'tag', new_tag, 'attribs', new_attribs, 'value', '', next, []);
      end
      str = str_t;
      eot = 0;
      %>
    else
      %< Appel du m�me code sur la suite

      % et stockage du resultat dans mes children.
      % Le code met aussi � jour le string courant |str|,
      % il en enl�ve la partie correspondant au string que je viens de trouver.
      [t,str] = parse_xml(str_t, new_tag, '', new_attribs, 1+idx);
      if isfield(t, next)
	nx = getfield( t, next);
        %nx = t.(next);
      else
	nx = [];
      end
      if isfield(z, next)
        nxt = getfield(z, next);
        nxt(end+1) = struct( 'tag', t.tag, 'attribs', t.attribs, 'value', t.value, next, nx);
        z   = setfield(z, next, nxt);
	%z.(next)(end+1) = struct( 'tag', t.tag, 'attribs', t.attribs, 'value', t.value, next, nx);
      else
	z = setfield(z, next, struct( 'tag', t.tag, 'attribs', t.attribs, 'value', t.value, next, nx));
	%z.(next) = struct( 'tag', t.tag, 'attribs', t.attribs, 'value', t.value, next, nx);
      end

      %>
    end
  end
  %>
end
%>

%< Parse attribs
function z =  parse_attribs( a)
if isempty(a)
  z = struct( 'name', '', 'value', '');
  return
end
b = tokens(a, ' ');
j = 1;
   
for i=1:length(b)
  if ~isempty(b{i})
    t = tokens(b{i}, '=');
    if length(t)==2
      u = t{2};
      if u(1)=='"'
	u = u(2:end);
      end
      if u(end)=='"'
	u = u(1:end-1);
      end
      z(j) = struct( 'name', (t{1}), 'value', u);
    else
        %here is the possible problem: if there is space between name and
        %'=' and value this will be counted 3 times as different names and
        %empty value
      z(j) = struct( 'name', (a), 'value', '');
    end
    j = j +1;
  end
end
%>


%<* Ecriture d'une structure xml
function z = write_xml(fid, xml_struct, idx)

next = 'child';

if nargin < 3
  idx = 0;
end

margin = repmat(' ',idx,1);

closed_tag = 1;
%< Ouverture du tag
if isfield(xml_struct, 'tag')
  closed_tag = 0;
  fprintf(fid, '%s<%s', margin, xml_struct.tag);
  %< Ecriture des attributs
  if ~isfield(xml_struct, 'attribs')
    error('xmltools:write_xml', 'malformed MATLAB xml structure : tag without attribs');
  end
  for i=1:length(xml_struct.attribs)
    if ~isempty(xml_struct.attribs(i).name)
        if xml_struct.tag(1) == '?'
            fprintf(fid, ' %s=%s', xml_struct.attribs(i).name, xml_struct.attribs(i).value);
        else
          fprintf(fid, ' %s="%s"', xml_struct.attribs(i).name, xml_struct.attribs(i).value);
        end
    end
  end
  %>
  
  %< Gestion des Auto closed tags
  % Si le tag n'est pas auto ferm�, alors |closed_tag| est � z�ro
  if ~isfield(xml_struct, next)
    error('xmltools:write_xml', 'malformed MATLAB xml structure : tag without %s', next);
  end
  if ~isfield(xml_struct, 'value')
    error('xmltools:write_xml', 'malformed MATLAB xml structure : tag without value');
  end
  if xml_struct.tag(1) == '?'
    fprintf(fid, '?>\n');
    closed_tag = 1;
  elseif isempty(getfield(xml_struct, next)) & isempty(xml_struct.value)
  %elseif isempty(xml_struct.(next)) & isempty(xml_struct.value)
    fprintf(fid, '/>\n');
    closed_tag = 1;
  elseif ~isempty(xml_struct.value)
      fprintf(fid, '>');
  else
      fprintf(fid, '>\n');
  end
  %>
end
%>

%< Ecriture de la value
if isfield(xml_struct, 'value')
  if ~isempty(xml_struct.value)
    fprintf(fid, '%s',  xml_struct.value);
  end
end
%>

%< Ecriture des enfants
if ~isfield(xml_struct, next)
  error('xmltools:write_xml', 'malformed MATLAB xml structure : tag without %s', next);
end
those_children = getfield(xml_struct, next);
%those_children = xml_struct.(next);

for i=1:length(those_children)
  write_xml(fid, those_children(i), idx+1);
end
%>

%< Fermeture du tag
if ~closed_tag
    if ~isempty(xml_struct.value)
        fprintf(fid, '</%s>\n',  xml_struct.tag);
    else
      fprintf(fid, '%s</%s>\n',  margin,xml_struct.tag);
    end
end
%>
%>*


%<* get childs with a specific tag name
function z = get_childs(z, next, tag_name);
u = getfield(z, next);
zo = [];
for i=1:length(u)
  v = u(i);
  if strcmp((v.tag), (tag_name))
    if isempty(zo)
      zo.anext= v;
    else
      zo.anext(end+1) = v;
    end
  end
end
if ~isstruct( zo)
  if isfield(z, 'tag')
    tn = z.tag;
  else
    tn = 'root?';
  end
  error('XMLTOOLS:GET-TEG', 'problem in finding tag <%s> under one <%s>', tag_name, tn);
end
z = [ zo.anext ];
%>*


%< udeblank
function s = udeblank(str)
s = deblank(str(end:-1:1));
s = s(end:-1:1);
if length(s)==0
  s = '';
end
%>

%< emptystruct
function z = emptystruct(next)
z = struct( 'tag', [], 'value', [], 'attribs', [], next, []);
%>

%< Tokens
function l = tokens(str,del)
l={} ; 
% Boucle sur les tokens.
del = sprintf(del) ;

space_pos = strfind(str,del);
eq_pos = strfind(str,'=');
index2squash = [];
for k=1:length(eq_pos)
    if sum(ismember(eq_pos(k)+[-1 1], space_pos))==2
        index2squash = [index2squash eq_pos(k)+[-1 1]];
    end
end
str(index2squash)='';

while ~isempty(str)
  [tok,str] = strtok(str,del) ;
  l{end+1} = tok ;
end

%fixing possible spaces btw '=' and name/value
% if length(l)==3 & strcmp(l{2},'=')
%     keyboard
% end

function [xml, rxml] = LoadXml(fbasename,varargin)

%function [xml, rxml] = LoadXml(FileBase)
%
% loads the xml file using xmltools (have to have it in the path)
% rxml returns it's original layout - very messy structure but contains all
% the xml file contents.
% xml - is the ouput structure which is backwards compatible to LoadPar
% output, so you can use it instead ..also loads some usefull stuff -
% Anatomoical groups with Skips , Spike electrode groups
% more can be added later (e.g. parameters of the process scripts)
% this script is written for xml version 1.1 .. older version doesn't work.
% additions are welcome

xml = struct;

xmli = strfind(fbasename,'.xml');
if isempty(xmli)
   fbasename = [fbasename '.xml'];
end
rxml = xmltools(fbasename);
try
    rxml = rxml.child(2);
catch
    rxml = rxml.child;
end
% from this level all children are the different parameters fields
xml.FileName = fbasename;

for i=1:length(rxml.child)

    switch lower(rxml.child(i).tag)
        
        case 'generalinfo'
            xml.Date = rxml.child(i).child(1).value; % date of xml file creation?

        case 'acquisitionsystem'
            xml.nBits = str2num(rxml.child(i).child(1).value); % number of bits of the file
            xml.nChannels = str2num(rxml.child(i).child(2).value);
            xml.SampleRate = str2num(rxml.child(i).child(3).value);
            xml.SampleTime = 1e6/xml.SampleRate; %to make backwards compatible
            xml.VoltageRange = str2num(rxml.child(i).child(4).value);
            xml.Amplification = str2num(rxml.child(i).child(5).value);
            xml.Offset =  str2num(rxml.child(i).child(6).value);
            
        case 'fieldpotentials'
            xml.lfpSampleRate = str2num(rxml.child(i).child.value);
            
        case 'anatomicaldescription'
            tmp = rxml.child(i).child.child;
            for grpI =1:length(tmp)
                for chI=1:length(tmp(grpI).child)
                    xml.AnatGrps(grpI).Channels(chI) = str2num(tmp(grpI).child(chI).value);
                    xml.AnatGrps(grpI).Skip(chI) = str2num(tmp(grpI).child(chI).attribs.value);
                end
            end
            
        case 'spikedetection'
            if ~isempty(rxml.child(i).child)
                tmp =rxml.child(i).child.child;
                for grpI =1:length(tmp)
                    for chI=1:length(tmp(grpI).child(1).child)
                        xml.SpkGrps(grpI).Channels(chI) = str2num(tmp(grpI).child(1).child(chI).value);
                    end
                    if length(tmp(grpI).child)>1
                        xml.SpkGrps(grpI).nSamples = str2num(tmp(grpI).child(2).value);
                        xml.SpkGrps(grpI).PeakSample = str2num(tmp(grpI).child(3).value);
                        xml.SpkGrps(grpI).nFeatures = str2num(tmp(grpI).child(4).value);
                    end
                    %backwards compatibility
                    xml.nElecGps = length(tmp);
                    xml.ElecGp{grpI} = xml.SpkGrps(grpI).Channels;
                end
            else
                xml.nElecGps = 0;
            end


        case 'programs'
            tmp = rxml.child(i).child;
            for i=1:length(tmp)
                if strcmp(tmp(i).child(1).value,'process_mhipass')
                    for j=1:length(tmp(i).child(2).child )
                        if strcmp(tmp(i).child(2).child(j).child(1).value,'frequency')
                            xml.HiPassFreq = str2num(tmp(i).child(2).child(j).child(2).value);
                            break
                        end
                    end
                end
            end
    end


end

function derivative = bz_Diff(samples,timestamps,varargin)

%Diff - Differentiate.
%
%  USAGE
%
%    derivative = Diff(samples,<options>)
%
%    samples        data to differentiate
%                   measured in number of samples (default = no smoothing)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'type'        'linear' if samples are linear values (default),
%                   'circular' otherwise
%     'smooth'      standard deviation for Gaussian kernel (default = 0)
%    =========================================================================
%
%  SEE
%
%    See also Interpolate.
%

% Copyright (C) 2004-2011 by Michaël Zugaro
% edited by David Tingley to fit buzcode formatting standards, 2017
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
type = 'linear';
smooth = 0;

if nargin < 1,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Diff">Diff</a>'' for details).');
end

if nargin < 1 | mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help Diff">Diff</a>'' for details).');
end

% Parse parameter list
for j = 1:2:length(varargin),
	if ~ischar(varargin{j}),
		error(['Parameter ' num2str(j+7) ' is not a property (type ''help <a href="matlab:help Diff">Diff</a>'' for details).']);
	end
	switch(lower(varargin{j})),
		case 'type'
			type = lower(varargin{j+1});
			if ~isastring(type,'linear','circular')
				error('Incorrect value for ''type'' (type ''help <a href="matlab:help Diff">Diff</a>'' for details).');
			end
		case 'smooth'
			smooth = varargin{j+1};
			if ~isdvector(smooth,'>=0') || length(smooth) > 2
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help Diff">Diff</a>'' for details).');
			end
        otherwise
			error(['Unknown property ''' num2str(varargin{j}) ''' (type ''help <a href="matlab:help Diff">Diff</a>'' for details).']);
	end
end

% Circular?
if strcmp(type,'circular')
	samples = unwrap(samples);
end

% Smooth
if smooth ~= 0
    samples = medfilt1(samples,smooth);
% 	if length(samples(1,2:end)) >= 2
% 		smooth = [smooth 0];
% 	end
% 	samples = Smooth(samples,smooth);
end

% Differentiate
derivative = zeros(size(samples));
% derivative(:,1) = samples(:,1);
t = timestamps;
dt = diff(t);
t_diff = t(1:end-1)+dt/2;
for i = 1:size(samples,2)
	d0 = diff(samples(:,i))./dt;
	d1 = interp1(t_diff,d0,t(2:end-1,1));
	derivative(:,i) = [d0(1);d1;d0(end)];
end
