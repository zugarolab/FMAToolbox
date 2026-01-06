function [cell_metrics,basepath] = ExportCellExplorer(basepath,varargin)

% brainRegions = session.brainRegions;

[~,basename] = fileparts(basepath);

filename = fullfile(basepath,[basename '.cell_metrics.cellinfo.mat']);
if exist(filename,'file')
    load(filename,'cell_metrics');
    return
end

[~,animalName] = fileparts(fileparts(basepath));

disp([datestr(clock) ': Starting session ' basepath]);
% for i=1:size(basepaths)
cd(basepath);
session = sessionTemplate(pwd,'showGUI',false,'basename',basename);
session.animal.name = animalName;
session.animal.sex = 'Male';
session.animal.species = 'Rat';
session.animal.strain = 'Long Evans';
session.brainRegions = brainRegions;

subsessions = LoadEvents(fullfile(basepath,[basename '.cat.evt']));
times = reshape(subsessions.time(:)',2,[])';
names = cellfun(@(x) x(14:end),subsessions.description(1:2:end),'UniformOutput',0); % After "beginning of"
for i=1:size(times,1)
    epoch = struct;
    epoch.name = names{i};
    epoch.startTime = times(i,1);
    epoch.stopTime = times(i,2);
    session.epochs{i} = epoch;
end

% Make channel coords
g = (Group(session.extracellular.electrodeGroups.channels{:}));
missing = find(~ismember((1:session.extracellular.nChannels)',g(:,1)));
session.extracellular.electrodeGroups.channels{end+1} = missing';
g = sortrows(Group(session.extracellular.electrodeGroups.channels{:}));
session.extracellular.chanCoords.x = g(:,2)*100;

y = sortby(sortby(CumSum(ones(size(g(:,1))),[0;diff(g(:,2))~=0]),g(:,1))*10,g(:,1));
g = sortrows(g);

try
    regions = []; fields = fieldnames(session.brainRegions);
    for i=1:length(fields)
        this = session.brainRegions.(fields{i})(:);
        this(:,2) = i;
        regions = [regions; this];
    end
    regions(end+1,:) = [size(regions,1)+1 0]; % add trash region
    x = g(:,2)+regions(g(:,2),2)*10; x = x*100;
catch
    x = g(:,2); x = x*100;
end

session.extracellular.chanCoords.x = x;
session.extracellular.chanCoords.y = -y;
session.extracellular.chanCoords.layout = 'hexatrodes';
save(fullfile(basepath,[basename '.session.mat']),'session');

spikes = loadSpikes('format','klustakwik','saveMat',true,'basepath',basepath,'basename',basename);
cell_metrics = ProcessCellMetrics('spikes',spikes,'session',session,'sessionSummaryFigure',false,'manualAdjustMonoSyn',false);

if isempty(cell_metrics.brainRegion{1})
    cell_metrics.brainRegion = fields(regions(cell_metrics.electrodeGroup,2))';
    saveAsFullfile = fullfile(basepath,[basenameFromBasepath(basepath),'.','cell_metrics','.cellinfo.','mat']);
    saveCellMetrics(cell_metrics,saveAsFullfile);
end
