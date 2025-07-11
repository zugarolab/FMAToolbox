function varargout = anovabar(data,groups,varargin)

%anovabar - Plots bars of the mean +/- s.e.m. for each group
%
%   Plots bars of the mean +/- s.e.m. for each group, indicating significant
%   differences between the groups. Each group is also tested against 0.
%   For non-parametric data, median +/- standard error of the median are 
%   displayed instead, and non-parametric tests are performed.
%
%  USAGE
%
%    [handles,h] = anovabar(data,groups,<options>)
%
%    data           data to test, provided either in matrix (each column is
%                   a different group) or a vector (grouping indicated through
%                   a grouping variable; see examples below)
%    group          grouping variable, indicating the group that each of the
%                   rows of 'data' corresponds to. Alternatively, it can be
%                   empty (default; only possible if 'data' is a matrix) or
%                   a string 'group' (indicating that the grouping variable
%                   is the last column of 'data')
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%    'alpha'        two values for the confidence levels for the following 
%                   tests: (1) between the data and 0, and (2) between
%                   groups (default = [0.05 0.05]). Set alpha to zero to 
%                   skip respective test (e.g. alpha = [0 0.05] would skip
%                   tests between each column and zero).
%    'parametric'   either 'on' (to display means and perform t-tests) 
%                   or 'off' (to display median and perform non-parametric 
%                   Wilcoxon signed-rank and ranksum tests / friedman).
%    'tail'         tail ('left', 'right', or 'both') for comparing each
%                   group to zero (default = 'both');
%    'correction'   type of statistical correction to use to control for
%                   multiple comparisons, i.e. between all pairs of multiple
%                   groups (default = tukey-kramer)
%    'precedence'   in case of a two-way anova, indicate which grouping variable
%                   should be tested while the other variable remains stable
%                   By default (default = 1), the first variable, i.e. the columns 
%                   of the matrix 'data' will be tested for a difference in 
%                   each of the groups described by the 2nd grouping variable.
%                   Set to 2 to inverse this behavior.                 
%    'paired'       a setting indicating whether values in the same row are
%                   paired (default = 'on').
%    'color'        bar colors, a matrix with three columns (RGB values) and
%                   where each line is a group
%    'baralpha'     transparency value for bars (between 0 and 1)
%    'barwidth'     the width of the bars to be plotted (default = 0.8).
%    'ns'           display the "n.s." label when between-group comparisons
%                   are not significant (default = false).
%    'labels'       x labels for conditions
%    'legend'       if 'off', plotted elements won't appear in legend (default =
%                   'on'); if different from 'on', must be a string array with
%                   group names
%    =========================================================================
%
%  OUTPUT
%
%    handles        a structure of handles for the objects in the figure.
%                   Includes the fields 'bar','errorbars','stars' (for the 
%                   stars above each bar indicating that values in that group
%                   are significantly different from 0), and 'comparisons'
%                   for all the objects (lines and stars) illustrating the 
%                   between-group comparisons.
%    h              h > 0 if the null hypothesis (that the groups are not 
%                   different) can be rejected (h values of 1, 2, and 3, 
%                   correspond to confidence levels 0.05, 0.01, and 0.001, 
%                   respectively)
%
%   Provide data as a matrix (each column will be treated separately), or
%   as grouped data with "group" as a grouping vector (each line indicating
%   which group the respective row in "data" corresponds to). Alternatively,
%   the grouping vector could be the last column of 'data' itself, in which case
%   you can call anovabar(data,'grouped');
%
%  EXAMPLES
%
%  Case 1: the two columns of 'data' correspond to two paired conditions
%  (e.g. control firing 2s before the intervals of interest and
%  the response of interest). Each row corresponds to paired observations:
%       data = [CountInIntervals(spikes(:,1),intervals-2) CountInIntervals(spikes(:,1),intervals)];
%       anovabar(data); % or anovabar(data,[]);
%  Note that for parametric data, the anova will not take the pairing into 
%  account (but for non-parametric data, a paired friedman test is performed).
%
%  Case 2: the observations are not paired; groups are therefore indicated
%  separately by a grouping variable:
%       controlData = CountInIntervals(spikes(:,1),baseline);
%       % note that here, 'baseline' intervals are not necessarily paired to the intervals of interest below
%       responseData = CountInIntervals(spikes(:,1),intervals);
%       data = [controlData; responseData]; groups = [ones(size(controlData)); ones(size(responseData))*2];
%       anovabar(data,groups);
%
%  Case 2bis: same as Case 2, but the grouping variable is contained within 'data':
%       controlData = CountInIntervals(spikes(:,1),baseline); controlData(:,end+1) = 1;
%       responseData = CountInIntervals(spikes(:,1),intervals); responseData(:,end+1) = 2;
%       data = [controlData; responseData];
%       anovabar(data,'grouped'); % or, of course, anovabar(data(:,1),data(:,end)) as in case 2
%
%  Case 3: two-way anova: data are grouped according to the columns in 'data'
%       AND ALSO by the grouping variable. In this case, the difference between the
%       (paired) columns will be tested for each of the groups indicated by the grouping
%       variable. For the opposite behavior (testing for the difference between the
%       (unpaired) groups indicated by the grouping variable, for each of the columns
%       of 'data', use <a href="matlab:help anovabar2">anovabar2</a>.
%
%   Inspired by barwitherr by Martina F. Callaghan
%
% Copyright (C) 2018-2022 by Ralitsa Todorova, 2025 by Pietro Bozzo (graphics)
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
alpha = [0.05 0.05];
correction = 'tukey-kramer'; 
grouped = false;
data = double(data);
parametric = true;
precedence = 1;
paired = true;
color = [];
barAlpha = 1;
barWidth = 0.8;
tail = 'both'; % tail for "test0" where each group is compared to zero
ns = false;
xlabels = string.empty;
legendLabels = "on";

for i = 1:2:length(varargin)
    if ~ischar(varargin{i})
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).']);
    end
    switch lower(varargin{i})
        case 'alpha'
            alpha = varargin{i+1};
            if ~isdvector(alpha) || length(alpha)>2
                error('Incorrect value for property ''alpha'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        case 'tail'
            tail = varargin{i+1};
            if ~isastring(lower(tail),'left','right','both')
                error('Incorrect value for property ''tail'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        case 'parametric'
            parametric = varargin{i+1};
            if isastring(lower(parametric),'on','off')
                if strcmpi(parametric,'on'); parametric = true; else, parametric = false; end
            end
            if ~islogical(parametric) || length(parametric)>1
                error('Incorrect value for property ''parametric'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        case 'paired'
            paired = varargin{i+1};
            if isastring(lower(paired),'on','off')
                if strcmpi(paired,'on'); paired = true; else, paired = false; end
            end
            if ~islogical(paired) || length(paired)>1
                error('Incorrect value for property ''paired'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        case 'correction'
            correction = varargin{i+1};
            if ~isastring(correction)
                error('Incorrect value for property ''correction'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        case 'precedence'
            precedence = varargin{i+1};
            if ~isdvector(precedence,'#1')
                error('Incorrect value for property ''precedence'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        case 'color'
            try % convert colors like 'k' to RGB triplet
                color = validatecolor(varargin{i+1},'multiple');
            catch ME
                throw(ME)
            end
        case 'baralpha'
            barAlpha = varargin{i+1};
            if ~isdvector(barAlpha,'#1','>=0','<=1')
                error('Incorrect value for property ''barAlpha'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        case 'barwidth'
            barWidth = varargin{i+1};
            if ~isdvector(barWidth,'#1')
                error('Incorrect value for property ''barWidth'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        case 'ns'
            ns = varargin{i+1};
            if isastring(lower(ns),'on','off')
                if strcmpi(ns,'on'); ns = true; else, ns = false; end
            end
            if ~islogical(ns) || length(ns)>1
                error('Incorrect value for property ''ns'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        case 'labels'
            xlabels = varargin{i+1};
            if ~isastring(xlabels)
                error('Incorrect value for property ''labels'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        case 'legend'
            legendLabels = string(varargin{i+1});
            if ~isvector(legendLabels)
                error('Incorrect value for property ''legend'' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).');
            end
        otherwise
            error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help anovabar">anovabar</a>'' for details).']);
    end
end

% if a single value for alpha is provided, use that value for both kinds of comparisons (to zero and between groups):
if length(alpha)==1, alpha = [alpha alpha]; end
if parametric
    test0 = @(x) out2(@ttest,x,[],'tail',tail);
    average = @(x) nanmean(x,1);
    semfun = @(x) nansem(x,1);
    testbetween = @anova1;
    if paired, testpaired = @(x) anova2(x,1,'off'); else testpaired = @(x) testbetween(x,[],'off'); end
else
    test0 = @(x) signrank(x,[],'tail',tail);
    average = @(x) nanmedian(x,1);
    semfun = @(x) semedian(x);
    testbetween = @kruskalwallis;
    if paired
        if size(data,2)==2, testpaired = @helper_signrank; 
        else, testpaired = @(x) friedman(x,1,'off'); 
        end
    else
        testpaired = @(x) testbetween(x,[],'off');
    end
end

if nargin==1 || isempty(groups)
    if sum(sum(~isnan(data))==0) > 0 % if one of the groups contains no real values
        return
    end
    xOrder = 1:size(data,2);
    values = average(data);
    errors = semfun(data);
else
    grouped = true;
    if isastring(groups,'groups','grouped','group')
        groups = data(:,end);
        data = data(:,1:end-1);
    end
    ok = ~isnan(groups);
    data = data(ok,:);
    groups = groups(ok);
    xOrder = unique(groups)';
    n = numel(xOrder);
    if precedence==1 % separate the groups provided by the grouping variable; test difference between columns of 'data' in the same group
        values = nan(n,size(data,2));
        errors = nan(n,size(data,2));
        for i=1:n
            group = xOrder(i);
            values(i,:) = average(data(groups==group,1:end));
            errors(i,:) = semfun(data(groups==group,1:end))';
        end
        if size(data,2)>1
            xOrder = 1:n;
        end
    else % separate the columns of 'data'; test difference between groups provided by the grouping variable in the same column
        values = nan(size(data,2),n);
        errors = nan(size(data,2),n);
        for i=1:n
            group = xOrder(i);
            values(:,i) = average(data(groups==group,1:end));
    		errors(:,i) = semfun(data(groups==group,1:end))';
    	end
    	if size(data,2)>1
    		xOrder = 1:size(data,2);
        end
    end
    groups = double(groups);
end

hStars0 = matlab.graphics.primitive.Line.empty;
hComp = matlab.graphics.chart.primitive.Line.empty;
hStars = matlab.graphics.primitive.Text.empty;

[nRows,nCols] = size(values);
% validate number of legend labels
if ~isscalar(legendLabels) && numel(legendLabels) ~= nCols
    error('Labels in ''legend'' must be as many as there are groups.')
end
% plot bars using built-in function
hbar = bar(xOrder,values,'BarWidth',barWidth,'FaceAlpha',barAlpha);
if ~isempty(color)
    for i = 1 : numel(hbar)
        hbar(i).FaceColor = color(i,:);
    end
end
hold on

vers = version;
if num2str(vers(1))<8 % Code for Matlab 2010
    if nRows > 1
        hErrorbar = zeros(1,nCols);
        for col = 1:nCols
            % Extract the x location data needed for the errorbar plots:
            x = get(get(hbar(col),'children'),'xdata');
            % Use the mean x values to call the standard errorbar function; the
            % errorbars will now be centred on each bar; these are in ascending
            % order so use xOrder to ensure y values and errors are too:
            hErrorbar(col) = errorbar(mean(x,1), values(xOrder,col), errors(xOrder,col), errors(xOrder, col), '.k');
            set(hErrorbar(col), 'marker', 'none')
        end
    else
        x = get(get(hbar,'children'),'xdata');
        hErrorbar = errorbar(mean(x,1), values, errors, errors, '.k');
        set(hErrorbar, 'marker', 'none');
    end
elseif verLessThan('matlab', '9.7') % Code for Matlab 2016 to 2019a, should we drop support for these few versions and just have > / < 2019b?
    xBars = vertcat(hbar.XData) + vertcat(hbar.XOffset); % xBars(i,j): ascissa of bar i in group j
    hErrorbar = errorbar(xBars.',values,errors,'Color','k','LineStyle','none');
else % Code for Matlab 2019b and above
    hErrorbar = errorbar(vertcat(hbar.XEndPoints).',values,errors,'Color','k','LineStyle','none');
end

%% Adding significance indication (stars etc)

% one way anova
if ~grouped || size(data,2)==1 
    if grouped, [~,~,stats] = testbetween(data(:,1),groups,'off'); u = unique(groups)';
    else [~,pTest,stats] = testpaired(data); u = 1:size(data,2);
    end
    if alpha(1)>0 % Unless alpha(1)=0, for each group, check if it's different from zero
        for i=1:length(u)
            if grouped, thesedata = data(groups==u(i),1); else thesedata = data(:,i); end
            if sum(~isnan(thesedata))>2
                p = test0(thesedata);
                if p<alpha(1)
                    % Put a little star above it to show it's significantly different from zero
                    hStars0(end+1,1) = plot(xOrder(i), sign(values(i))*(max([abs(values(i)+errors(i)) abs(values(i)-errors(i))])+abs(values(i)/7)), '*', 'markersize', 5, 'MarkerEdgeColor', [0 0 0.5]);
                end
            end
        end
    end
    u = unique(xOrder)';
    sigFor2Groups = nan(length(u),1);
    if alpha(2)>0 % Unless alpha(2)=0, perform tests between groups
        if ischar(stats) && strcmp(stats,'signrank') % signrank test doesn't give a "stats" output compatible with multcompare
            comparison = [1 2 (pTest<alpha(2)) + (pTest<alpha(2)/5) (pTest<alpha(2)/50)];
        else
            comparison = multcompare(stats,'display', 'off','alpha',alpha(2),'ctype',correction); % default alpha = 0.05
            comparison = [comparison(:,1:2) double(comparison(:,3).*comparison(:,5)>0)]; %if the upper and lower bound have the same sign
            comparison2 = multcompare(stats,'display', 'off', 'alpha', alpha(2)/5,'ctype',correction); % **, default alpha = 0.01
            comparison(:,3) = comparison(:,3) + double(comparison2(:,3).*comparison2(:,5)>0); % third column shows the number of stars to be included. 1 for 0.05, 1 more for 0.01, and another one for 0.001       if sum(double(comparison2(:,3).*comparison2(:,5)>0)),
            comparison3 = multcompare(stats,'display', 'off', 'alpha', alpha(2)/50,'ctype',correction); % ***, default alpha = 0.001
            comparison(:,3) = comparison(:,3) + double(comparison3(:,3).*comparison3(:,5)>0);
        end
        % Now to plot them:
        for i=1:size(comparison,1)
            s = sign(Portion(values(:)>0)./(Portion(values(:)>0 | values(:)<0))-0.5); if s==0,s=sign(nanmean(values(:)));end
            smallnumber = nanmean(values(:))/10; %for display purposes, so things don't overlap
            y1 = s*max(abs([errors(:)+values(:); -errors(:)+values(:)]))+(comparison(i,2)-comparison(i,1))*smallnumber;
            y2 = y1+smallnumber/4;
            x1 = xOrder(comparison(i,1)) - (xOrder(comparison(i,2))-xOrder(comparison(i,1)))*0.2 + 0.4;
            x2 = xOrder(comparison(i,2)) + (xOrder(comparison(i,2))-xOrder(comparison(i,1)))*0.2 - 0.4; %lines between adjecent bars are shorter; the bigger the line, the farther it goes.
            if comparison(i,3)==0
                if ns, text(mean([comparison(i,1) comparison(i,2)]), y2, 'n.s.'); end
            else
                hComp(end+1,1) = plot([x1 x2], [y1 y1], 'k');
                hComp(end,2) = plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                hComp(end,3) = plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                hStars(end+1,1) = text(mean([xOrder(comparison(i,1)) xOrder(comparison(i,2))]), y2, repmat('*',1,comparison(i,3)));
            end
        end
    end
    comparison = sigFor2Groups;

    % two way anova
else
    if num2str(vers(1))<8 % Code for Matlab 2010
        x = get(hbar,'children');
        for i=1:length(x)
            xs(:,i) = nanmean(get(x{i},'xdata'));
        end
    else % New code for Matlab 2016
        for i=1:length(hbar)
            xs(:,i) = hbar(i).XData + hbar(i).XOffset;
        end
    end
    u = unique(groups)';
    % Test if each group is different from zero and draw stars
    if alpha(1)>0 % If alpha(1)=0, skip this
        for j=1:length(u), for i=1:size(data,2)
                if precedence==1, thisx = xs(j,i); thisy = sign(values(j,i))*(max([abs(values(j,i)+errors(j,i)) abs(values(j,i)-errors(j,i))])+abs(values(j,i)/7));
                else thisx = xs(i,j); thisy = sign(values(i,j))*(max([abs(values(i,j)+errors(i,j)) abs(values(i,j)-errors(i,j))])+abs(values(i,j)/7));
                end
                thesedata = data(groups==u(j),i);
                if sum(~isnan(thesedata))>2, p = test0(thesedata);
                    if p<alpha(1)
                        hStars0(end+1,1) = plot(thisx, thisy, 'k*', 'markersize', 5, 'MarkerEdgeColor', [0 0 0]);
                    end % Put a little star above it
                end
            end; end
    end
    % Test if groups are different between each other
    sigFor2Groups = nan(length(u),1);
    if alpha(2)>0 % If alpha(2)=0, skip this
        if precedence==1, maxj = length(u); else, maxj = length(xOrder); end
        for j=1:maxj
            if precedence==1, [~,pTest,stats] = testpaired(data(groups==u(j),:));
            else, [~,~,stats] = testbetween(data(:,j), groups,'off');
            end
            if ischar(stats) && strcmp(stats,'signrank') % signrank test doesn't give a "stats" output compatible with multcompare
                comparison = [1 2 (pTest<alpha(2)) + (pTest<alpha(2)/5) (pTest<alpha(2)/50)];
            else
                try
              		 comparison = multcompare(stats,'display', 'off','alpha',alpha(2),'ctype',correction);
                     comparison = [comparison(:,1:2) double(comparison(:,3).*comparison(:,5)>0)]; %if the upper and lower bound have the same sign
                     comparison2 = multcompare(stats,'display', 'off', 'alpha', alpha(2)/5,'ctype',correction);
                     comparison(:,3) = comparison(:,3) + double(comparison2(:,3).*comparison2(:,5)>0); % third column shows the number of stars to be included. 1 for 0.05, 1 more for 0.01, and another one for 0.001       if sum(double(comparison2(:,3).*comparison2(:,5)>0)),
                     comparison3 = multcompare(stats,'display', 'off', 'alpha', alpha(2)/50,'ctype',correction);
                     comparison(:,3) = comparison(:,3) + double(comparison3(:,3).*comparison3(:,5)>0);
                catch
                    continue % if a whole column of data is missing, then comparisons are impossible
                end
            end
            sigFor2Groups(j,1) = comparison(1,3);
            for i=1:size(comparison,1)
                if precedence==1, s = sign(nanmean(nanmean(data(groups==u(j),:)))); else, s = sign(nanmean(data(:,j))); end
                smallnumber = nanmean(abs(values(:)))/10; %for display purposes, so things don't overlap
                y1 = s*max(abs([errors(j,:)'+values(j,:)'; -errors(j,:)'+values(j,:)']))+(comparison(i,2)-comparison(i,1))*smallnumber;
                y2 = y1+smallnumber/3;
                x1 = xs(j,comparison(i,1));
                x2 = xs(j,comparison(i,2));
                if comparison(i,3)==1
                    hComp(end+1,1) = plot([x1 x2], [y1 y1], 'k');
                    hComp(end,2) = plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                    hComp(end,3) = plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                    hStars(end+1,1) = text(mean([x1 x2]), y2, '*','HorizontalAlignment','center');
                elseif comparison(i,3)==2
                    hComp(end+1,1) = plot([x1 x2], [y1 y1], 'k');
                    hComp(end,2) = plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                    hComp(end,3) = plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                    hStars(end+1,1) = text(mean([x1 x2]), y2, '**','HorizontalAlignment','center');
                elseif comparison(i,3)==3
                    hComp(end+1,1) = plot([x1 x2], [y1 y1], 'k');
                    hComp(end,2) = plot([x1 x1], [y1 y1-smallnumber/10], 'k');
                    hComp(end,3) = plot([x2 x2], [y1 y1-smallnumber/10], 'k');
                    hStars(end+1,1) = text(mean([x1 x2]), y2, '***','HorizontalAlignment','center');
                end
            end
        end
    end
    comparison = sigFor2Groups;
end
hold off

% make x labels
if ~isempty(xlabels)
    xticks(hbar(1).XData)
    xticklabels(xlabels)
end

% handle legend
if ~isscalar(legendLabels) || legendLabels ~= "on"
    if legendLabels == "off"
        RemoveFromLegend(hbar)
    else
        for i = 1 : numel(hbar)
            hbar(i).DisplayName = legendLabels(i);
        end
    end
end
if ~isscalar(legendLabels) || legendLabels == "off"
    RemoveFromLegend(hErrorbar)
    RemoveFromLegend(hComp)
end

if nargout>0
    handles.bar = hbar;
    handles.errorbars = hErrorbar;
    handles.stars = hStars0;
    handles.comparisons = [hComp,hStars];
    varargout{1} = handles;
    varargout{2} = comparison;
end

% ------------------------------- Helper functions -------------------------------

function [h,p,stats] = helper_signrank(data,varargin)

ranks = data; ranks(:) = tiedrank(data(:));
p = signrank(diff(data,[],2));
h = p<0.05;
stats = 'signrank';

