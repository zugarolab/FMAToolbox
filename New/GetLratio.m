function [lratio, distance] = GetLratio(units)

%Lratio - Determine the isolation quality of one or more clusters.
%
% Compute isolation distance, 
%
%  USAGE
%
%    lratio = Lratio(units)
%
%    units          list of units to process (default = all clusters except 
%                   0 and 1)
%    
%  OUTPUT
% 
%    lratio          list of lratio for each cluster
%    distance        Distance obtain with IsolationDistance.m, as defined 
%                    by K. Harris (see Schmitzer-Torbert et al. 2005). : 
%                    For each cluster, a measure of the separation between 
%                    the cluster and all other spikes. Namely, the 
%                    Mahalanobis distance of the Nth closest spike outside 
%                    the cluster (where N is the number of spikes inside 
%                    the cluster).
%
% Lratio = (sum(i~=C)(1-CDFqui2df=4(D2i,c)))/nb spk on this electrod group;
%
% Céline Drieu
%
% Note: The L-ratio and isolation distance computed here take into account all other 
% spikes on the same group, including spikes in noise clusters.

values1 = nan(size(units,1),1);
values2 = nan(size(units,1),1);

groups = unique(units(:,1));

n=1;
% Loop through electrod groups
for i=1:length(groups),
    
    features = GetSpikeFeatures(groups(i));
    nFeatures = size(features,2)-3;
    
    unitsgroup = units(units(:,1)==groups(i),:);
    nUnitsgroup = size(unitsgroup,1);
    
    spikesInElectrodeGroup = features(:,2)==groups(i);
    nSpikesInElectrodeGroup = sum(spikesInElectrodeGroup);
    
    % Loop through units
    for j = 1:nUnitsgroup,
        unit = unitsgroup(j,:);     
        unitStr = ['(' int2str(unit(1)) ',' int2str(unit(2)) ')'];

        % Determine spikes inside vs outside cluster
        spikesInCluster = features(:,2)==unit(1) & features(:,3)==unit(2);
        nSpikesInCluster = sum(spikesInCluster);
        
        % We need at least as many spikes as features, but no more than half the total number of spikes
        if nSpikesInCluster < nFeatures,
            warning(['Not enough spikes for unit ' unitStr '.']);
            continue
        end
        if nSpikesInCluster > nSpikesInElectrodeGroup/2,
            warning(['Too many spikes for unit ' unitStr '.']);
            continue
        end
        % Compute Mahalanobis distance for spikes outside (mOut) the cluster
        m = mahal(features(:,4:end),features(spikesInCluster,4:end));
        mOut = m(~spikesInCluster);
        
        % Determine the Mahalanobis of the Nth closest spike outside the cluster (where N is the number of spikes inside the cluster)
        sOut = sort(mOut);
        distance(n) = sOut(nSpikesInCluster);
        
        % Cumulative distribution function of chi2 distribution
        df = 4; % degrees of freedom
        p1 = cdf('chi2',mOut,df);
        
        % Compute the Lratio (normalized the sum by the total numberspikes
        % recorded on the tetrode, like in Pfeiffer&Foster,2013)
        lratio(n) = sum(1-p1)/nSpikesInElectrodeGroup;
        n=n+1;
    end
end





