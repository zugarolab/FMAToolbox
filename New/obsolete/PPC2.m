function ppc = PPC2(phases1,phases2)

% Pairwise Phase Consistency as defined by  
% Vincka, Wingerden, Womelsdorf, Fries, Pennartz (2010, NeuroImage)

difference = bsxfun(@minus,phases1,phases2');
ppc = nanmean(cos(difference(:))); % as cos(a-b) = cos(a)*cos(b)+sin(a)*sin(b), Vinck et al's definition of PPC

% n1 = length(phases1);
% n2 = length(phases2);
% difference = bsxfun(@minus,phases1,phases2');
% d = abs(wrap(difference));
% D = sum(d(:))/(n1*n2);
% 


