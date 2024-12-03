function data = logfactorial(data)

% We compute log(n!) as the sum of logs, i.e. log(n!) = sum log(i) for i=1:n
% First determine the largest n in the array
m = max(data(:));
% Create a look-up vector of sum log(i) for each i up to the largest n
sums = [0 cumsum(log(1:m))];
% Look-up the value for each item in the array
data(:) = sums(data+1);
end
