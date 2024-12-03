function samples = Sample(histogram, n)

% sample from distributions.
% EXAMPLE
% [h,t] = hist(spikes(:,1),resolution);
% samplespikes = t(Sample(h,1000));
% NOTE: multiple histograms may be provided in different columns of h.
%
% By Raly Todorova (2016), inspired by Dahua Lin's discretesample (2008)

if nargin<2
    error('Please provide 2 arguments: pdf(s)/histogram(s) and n of desired samples');
end

if isvector(histogram)
    histogram = histogram(:);
end

nDistributions = size(histogram,2);
nBins = size(histogram,1);

histogram = bsxfun(@rdivide,histogram,nansum(histogram));

% The main reasoning of this code was inspired by Dahua Lin's discretesample (2008)
% The reasoning is as follows:
% We can sample from a uniform distribution using rand and then check where it falls relative to some intervals
% We have nBins intervals. If rand(i) falls within the 16th bin, sample(i)=16.
% Now, to adjust the probability of obtaining each possible value, we can adjust the interval durations to reflect out distribution:
% The i-th interval duration should be proportional to the pdf (histogram) value at the i-th bin provided.
% This is achieved by creating intervals that start at 0 and build up with the values of the pdf (histogram):
cdf = [zeros(1,nDistributions); cumsum(histogram)];
intervals = [reshape(cdf(1:end-1,:),[],1) reshape(cdf(2:end,:),[],1)];
% the intervals need to be sorted, or else CountInIntervals skips them
[intervals,order] = sortrows(intervals);
count = CountInIntervals(sort(rand(n,1)),intervals);
% get a translation vector that would reverse the ordering of the intervals
inverse(order,1) = 1:length(order);
count = count(inverse,1);

% The generalisation to multiple distributions is achieved by having them all stacked below one another counts = [countsFor1; countsFor2;...;countsForLast];
% Sampling is achieved by repeating the i-th id count(i) times (e.g. rand(i) falls within the 16th bin => count(16) is increased by 1, and 16 is sampled
samples = repelem(rem((1:nBins*nDistributions)-1,nBins)'+1,count);
% Reshape the final result so that the samples stacked below one another are ordered by respective distributions (one per column)
samples = reshape(samples,n,nDistributions);
samples = Shuffle(samples')';

end

%% HELPER FUNCTIONS

function shuffled = Shuffle(matrix,Logical)

% row by row
% if logical is provided, it should indicate which elements are to shuffle (0s will remain as given)

if nargin>1
    [~,indices] = sort(rand(size(matrix')).*(-Logical'));
else
    [~,indices] = sort(rand(size(matrix')));
end
shuffled = matrix(sub2ind(size(matrix),meshgrid(1:size(matrix,1),1:size(matrix,2))',indices'));

end
