function ppc = PPCtrials(phases,trial)

% P(hat)2 Pairwise Phase Consistency as defined by  
% Vinck, Battaglia, Womelsdorf, Pennartz (2012, J Comp Neurosci),

phases = phases(:); trial = (:);
if numel(trial)~=numel(phases),
    error('Please provide a trial number for each phase');
end

nTrials = max(trial);
ppc = nan(nTrials);
for m=1:nTrials,
    for l=1:nTrials,
        if l~=m,
            ppc(m,l) = PPC2(phases(trial==m),phases(trial==l));
        end
    end
end

% as the diagonal is now composed of nans, there are n*(n-1) non-nan elements in the matrix 'ppc' 
% i.e. taking the mean is equilalent to the sum, normalised by n*(n-1) as in Vinck et al (2012)
ppc = nanmean(ppc(:));