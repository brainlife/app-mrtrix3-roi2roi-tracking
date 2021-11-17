function [] = classificationGenerator()

if ~isdeployed
    disp('loading path')
    addpath(genpath('/N/u/hayashis/git/vistasoft'))
    addpath(genpath('/N/u/brlife/git/jsonlab'))
    addpath(genpath('/N/u/brlife/git/wma_tools'))
end

% Load configuration file
config = loadjson('config.json')
roiPair = strtrim(config.roiPair)

% Set tck file path/s
disp('merging tcks')
tcks=dir('track*.tck')
jj=1;
kk=1;
for ii = 1:length(tcks)
    fgtmp = fgRead(fullfile(tcks(ii).folder,tcks(ii).name));
    if length(fgtmp.fibers) > 0
        fgPath{jj} = tcks(ii).name;
        jj=jj+1;
    else
        display(sprintf('track %s has 0 streamlines',tcks(ii).name))
        missing_tracts{kk} = tcks(ii).name;
        kk=kk+1;
    end
end
disp(fgPath)
[mergedFG, classification]=bsc_mergeFGandClass(fgPath);
%fgWrite(mergedFG, 'track/track.tck', 'tck');

if ~exist('wmc', 'dir')
    mkdir('wmc')
end
if ~exist('wmc/tracts', 'dir')
    mkdir('wmc/tracts')
end

% Amend name of tract in classification structure
roiPair = split(roiPair);
% for ii = 1:length(roiPair)/2
%     classification.names{ii} = strcat('ROI_',roiPair{(2*ii) - 1},'_ROI_',roiPair{(2*ii)});
% end
for ii = 1:length(classification.names)
    jj = (2*str2num(extractAfter(classification.names{ii},'track')))-1;
    kk = str2num(extractAfter(classification.names{ii},'track'))*2;
    roi_names = [roiPair(jj),roiPair(kk)];
    classification.names{ii} = strcat('ROI_',roi_names{1},'_ROI_',roi_names{2});
end
save('wmc/classification.mat','classification')
save('wmc/missing_tracts.mat','missing_tracts')

% split up fg again to create tracts.json
fg_classified = bsc_makeFGsFromClassification_v4(classification,mergedFG);
tracts = fg2Array(fg_classified);
%cm = parula(length(tracts));
cm = distinguishable_colors(length(tracts));
for it = 1:length(tracts)
   tract.name   = strrep(tracts{it}.name, '_', ' ');
   all_tracts(it).name = strrep(tracts{it}.name, '_', ' ');
   all_tracts(it).color = cm(it,:);
   tract.color  = cm(it,:);

   %tract.coords = tracts(it).fibers;
   %pick randomly up to 1000 fibers (pick all if there are less than 1000)
   fiber_count = length(tracts{it}.fibers);
   tract.coords = tracts{it}.fibers; 
   
   savejson('', tract, fullfile('wmc','tracts', sprintf('%i.json',it)));
   all_tracts(it).filename = sprintf('%i.json',it);
   clear tract
end

% Save json outputs
savejson('', all_tracts, fullfile('wmc/tracts/tracts.json'));

% Create and write output_fibercounts.txt file
for ii = 1 : length(fg_classified)
    name = fg_classified{ii}.name;
    num_fibers = length(fg_classified{ii}.fibers);
    
    fibercounts(ii) = num_fibers;
    tract_info{ii,1} = name;
    tract_info{ii,2} = num_fibers;
end

T = cell2table(tract_info);
T.Properties.VariableNames = {'Tracts', 'FiberCount'};

writetable(T, fullfile('wmc','output_fibercounts.txt'));

exit;
end



