function [width, trough2peak] = GetWaveformStats(waveforms)

% Stats are in index values. To convert them to seconds, multiply with your sampling frequency

% Check that we have a single waveform
if isvector(waveforms),
    waveforms = waveforms(:)';
end

n = size(waveforms,1);
for i=1:n,
    try
        waveform = waveforms(i,:);
        waveform = waveform - waveform(1); % remove baseline constant, so half to trough is trough / 2
        
        [troughValue,trough] = min(waveform);
        [~,trough2peak(i,1)] = max(waveform(trough:end));
        
        % the spike width is the distance between the two points at half-trough value surrounding the trough itself
        start = find(waveform < troughValue/2,1,'first'); % first point
        width(i,1) = find(waveform(start:1:end) > troughValue/2,1,'first'); % distance to the end
        
    catch
        width(i,1)=nan;
        trough2peak(i,1)=nan;
    end
end