function bursts = FindBursts(spikes,varargin)

% This function implements the method used by Wierzynski 2009. I.e. smooths
% the mua spike train and identified periods where it surpasses a low-threshold,
% provided that a high threshold is surpassed somewhere within the burst.


dt = 0.05; % in s ; 50 ms in Wierzynski
threshold = 2; %in sd; 2 in in Wierzynski
lowThreshold = 2;
smooth = (0.05/3)/dt; % 3σ = 50 ms in Wierzynski

for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'threshold',
            threshold = varargin{i+1};
            if ~isdvector(threshold,'#1','>0'),
                error('Incorrect value for property ''threshold'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        case 'lowthreshold',
            lowThreshold = varargin{i+1};
            if ~isdvector(lowThreshold,'#1','>0'),
                error('Incorrect value for property ''lowThreshold'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        case 'window'
            dt = varargin{i+1};
            if ~isdvector(dt,'#1','>0'),
                error('Incorrect value for property ''window'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        case 'smooth',
            smooth = varargin{i+1}/dt;
            if ~isdvector(smooth,'#1','>0'),
                error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).');
            end
        otherwise,
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindBursts">FindBursts</a>'' for details).']);
    end
end

spikes = spikes(:,1);
[binned,t] = binspikes(spikes(:,1),1/dt); t = t(:);
binned = Smooth(binned,smooth);

z = zscore(binned);

% Find periods depassing the low threshold
bursts(:,[1 3]) = t(ToIntervals(z>lowThreshold));

% Find maximal z value within putative bursts
[in,id] = InIntervals(t,bursts(:,[1 3]));
[peak,idx] = Accumulate(id(in),z(in),'mode','max');
t = t(in);
bursts(:,2) = t(idx);

% Apply high threshold
pass = peak>threshold;
bursts = bursts(pass,:);

end


%% HELPER FUNCTIONS
% Chronux' binspikes
function [dN,t]=binspikes(data,Fs,t)
% bin spikes at a specified frequency sampling i.e. sampling rate 1/sampling
% eg: 1ms accuracy use sampling = 1000
% Usage: [dN,t]=binspikes(data,Fs,t)
% Inputs:
% data   (data as a structure array of spike times; or as a single
%        vector of spike times)
% Fs     (binning frequency)
% t      (the minimum and maximum times to be used to form the bins - [mint maxt]
%            - optional. Default use the spike times themselves to
%              determine the location of the bins.
% Note: the times in data can be in any units. However, it is important
% that all units are chosen consistently. So, if spike times are in secs,
% Fs and t (if present) have to be in Hz and secs respectively. If spike
% times are in number of samples, Fs has to be 1, and t has to be in number
% of samples.
% Outputs:
% dN     (output binned spike counts as a matrix defined on bins starting with the
%         earliest spike across all channels and ending with the latest spike)
% t      (lower limit of each bin)
if nargin < 2; error('Need at least two input arguments'); end;
dt=1/Fs;
dtmp='';
if isstruct(data);
    C=length(data);
    fnames=fieldnames(data);
    if nargin <3 || isempty(t);
        mintime=zeros(1,C);
        maxtime=zeros(1,C);
        for ch=1:C
            eval(['dtmp=data(ch).' fnames{1} ';'])
            mintime(ch)=min(dtmp);
            maxtime(ch)=max(dtmp);
        end
        mintime=min(mintime);
        maxtime=max(maxtime);
    else
        %        maxtimech=zeros(1,C);
        %        for ch=1:C
        %          eval(['dtmp=data(ch).' fnames{1} ';'])
        % %          mintimech(ch)=min(dtmp);
        %          maxtimech(ch)=max(dtmp);
        %        end
        mintime=t(1);
        maxtime=t(end);
        %        mintimech=min(mintimech);
        %        maxtimech=max(maxtimech);
        %        if maxtimech > max(t); t=[t maxtimech+dt]; end;
    end
    t=linspace(mintime,maxtime,1+(maxtime-mintime)/dt);
    for ch=1:C;
        eval(['dtmp=data(ch).' fnames{1} ';'])
        x=histc(dtmp,t);
        dN(:,ch)=x(:);
    end
else
    dtmp=data;
    if nargin < 3;
        mintime=min(dtmp);
        maxtime=max(dtmp);
    else
        mintime=t(1);
        maxtime=t(end);
    end
    t=linspace(mintime,maxtime,1+(maxtime-mintime)/dt);
    if max(dtmp)>max(t); t=[t maxtime+dt]; end;
    x=histc(dtmp,t);
    dN=x(:);
end
end