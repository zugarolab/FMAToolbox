function [cosine,p,z,count] = BiasCosine(ibeforej1,count1,ibeforej2,count2,nIterations)

% as defined by Zachary Roth https://arxiv.org/pdf/1603.02916.pdf, page 65-67
% Note that bias matrices should be flattened (first dimension should be flattened bias,
% second dimension is sequence ID).
% Note that we are comparing sequence2 to sequence1, so we only take the variability of sequence 2 into account!!

if ~exist('nIterations','var'), nIterations = 100; end
if nargin<3, ibeforej2 = ibeforej1; count2 = count1; end

% norms1 = sqrt(sum(bias1.^2));
% norms2 = sqrt(sum(bias2.^2));

%%
nUnits = size(ibeforej1,1);
bias1 = ibeforej1 ./ count1;
bias1 = reshape(bias1,nUnits^2,[]);

ibeforej2 = reshape(ibeforej2,nUnits^2,[]); count2 = reshape(count2,nUnits^2,[]);
bias2 = ibeforej2 ./ count2;

% allow for negative biases:
bias1 = bias1*2-1;
bias2 = bias2*2-1;

notnan1 = double(~isnan(bias1)); notnan2 = double(~isnan(bias2));
bias1(~notnan1) = 0; bias2(~notnan2) = 0;
norms = (sqrt((bias1.^2)'*notnan2).*sqrt(notnan1'*(bias2.^2)));

cosine = bias1'*bias2 ./ norms;
count = notnan1'*notnan2;

%%
if nargout<2
    return
end
[p,z] = deal(nan(size(cosine)));
if nIterations==0
    return
end
count = notnan1'*notnan2;
norm1 = (sqrt((bias1.^2)'*notnan2));
try
    if log10(numel(cosine))>7, error('this requires too much memory'); end
    tic; shuffled = nan([size(cosine) nIterations]);
    if toc>30, error('this requires too much memory'); end
    [i,j] = meshgrid(1:size(count2,2),1:size(count2,1)); i=i(:); j=j(:);
    i = repelem(i,count2(:)); j = repelem(j,count2(:));
    for iteration = 1:nIterations,
        r = rand(sum(count2(:)),1)>0.5; % order between any 2 spikes is random, a flip of a coin
        ibeforej2shuffled = Accumulate([j i],r,'size',size(count2));
        bias2shuffled = (ibeforej2shuffled ./ count2)*2-1;
        notnan2 = double(~isnan(bias2shuffled));
        bias2shuffled(~notnan2) = 0;
        norms = norm1.*sqrt(notnan1'*(bias2shuffled.^2));
        shuffled(:,:,iteration) = bias1'*bias2shuffled ./ norms;
        if rem(iteration,50)==0,
            %             display([iteration round(toc)]);
            tic; end
    end
    p = mean(bsxfun(@le,abs(cosine),abs(shuffled)),3);
    z = (cosine - mean(shuffled,3))./var(shuffled,[],3);
catch
    for i=1:size(bias2,2),tic;
        idx = repelem((1:nUnits^2)',count2(:,i));
        shuffled = nan([size(cosine,1) nIterations]);
        for iteration = 1:nIterations,
            r = rand(sum(count2(:,i)),1)>0.5; % order between any 2 spikes is random, a flip of a coin
            ibeforej2shuffled = accumarray(idx,r,[nUnits^2 1]);
            bias2shuffled = (ibeforej2shuffled ./ count2(:,i))*2-1;
            bias2shuffled(~notnan2(:,i)) = 0;
            norms = norm1(:,i).*sqrt(notnan1'*(bias2shuffled.^2));
            shuffled(:,iteration) = bias1'*bias2shuffled ./ norms;
        end
        p(:,i) = mean(bsxfun(@le,abs(cosine(:,i)),abs(shuffled)),2);
        z(:,i) = (cosine(:,i) - mean(shuffled,2))./var(shuffled,[],2);
        if rem(i,1000)==0,display([i size(bias2,2)]); toc;tic; end
    end
end

% === HELPER FUNCTION ===
function [matrix,ibeforej,count] = GetBiasMatrix(sequenceList,nUnits)

nSequences = size(sequenceList,1);
%% Get a matrix descriptive of the firing statistics of neurons (how often neuron i fires after neuron j)
ibeforej = zeros(nUnits,nUnits,nSequences);
count = zeros(nUnits,nUnits,nSequences);
for i=1:nUnits,
    for j=i+1:nUnits,
        ok = find(sum(sequenceList==i,2)>0 & sum(sequenceList==j,2)>0);
        [iRank,iID] = find(sequenceList(ok,:)'==i);
        [jRank,jID] = find(sequenceList(ok,:)'==j);
        instances = Accumulate(iID,1,length(ok)).*Accumulate(jID,1,length(ok));
        count(i,j,ok) = instances;
        once = instances==1;
        if sum(once)>0,ibeforej(i,j,ok(once)) = iRank(ismember(iID,find(once)))<jRank(ismember(jID,find(once)));end
        for k=find(~once)',ibeforej(i,j,ok(k)) = sum(sum(bsxfun(@lt,iRank(iID==k),jRank(jID==k)')));end
    end
end
jbeforei = count-ibeforej;

%fill in lower part as well
count = count+permute(count,[2 1 3]);
ibeforej = ibeforej + permute(jbeforei,[2 1 3]);

matrix = ibeforej./count;

