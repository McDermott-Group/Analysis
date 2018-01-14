function [data]=plotDirSequenceFidelities(loadedDirData)

% Select a directory if loadedDirData doesn't exist and load the data from
% each file.
if ~exist('loadedDirData', 'var')
    
    pathname = '';
    fid = fopen(fullfile(tempdir, 'plotMeasurementData_last_pathname.txt'),...
      'r');
    if fid ~= -1
        pathname = fgetl(fid);
        fclose(fid);
    end
    
    starting_dir = fullfile('Z:\','mcdermott-group','data','sfq','Wafer070717C','C2','DR1');
    if strcmp(pathname, '') || ~exist(pathname, 'dir')
        if exist(starting_dir, 'dir')
            pathname = starting_dir;
        end
    end
    dataDir = uigetdir(pathname,'Select a folder with calibrated blob sequences');
    
    dataFiles = dir(fullfile(dataDir,'*.txt'));
    for file=dataFiles'
        data.sets.(strrep(file.name,'.txt',''))=loadMeasurementData(fullfile(dataDir,file.name));
    end
    data.dataDir = dataDir;
    
    % Save the last path used.
    fid = fopen(fullfile(tempdir,...
        'plotMeasurementData_last_pathname.txt'), 'w');
    if fid ~= -1
        fprintf(fid, '%s', pathname);
        fclose(fid);
    end
else
    data = loadedDirData;
end

% With the data loaded, grab the occupations and gate sequences of each run
data.seq_fids = zeros(length(data),1);
data.gate_seqs{length(data),1} = [];
filenames = fieldnames(data.sets);

for k=1:length(filenames)
   data.seq_fids(k) = 1 - data.sets.(filenames{k}).Poor_Mans_Occupation;
   data.gate_seqs{k} = data.sets.(filenames{k}).SFQ_User_Defined_Sequence;
end

% now plot it with its mean and std on the same plot
createFigure; hold on;
h(1)=bar(data.seq_fids,'c');
data.seq_fid_mean = mean(data.seq_fids);
data.seq_fid_std = std(data.seq_fids);
h(2)=refline(0,data.seq_fid_mean);
h(2).Color = 'k';
hline=refline(0,data.seq_fid_mean+data.seq_fid_std);
hline.Color = 'r';
hline=refline(0,data.seq_fid_mean-data.seq_fid_std);
hline.Color = 'r';
ax=gca;
set(ax,'XTickLabelRotation',300)
set(ax,'XTick',1:1:length(filenames))
set(ax,'XTickLabel',data.gate_seqs)
set(ax,'YGrid','on')
xlabel('Gate Sequence')
ylabel('Sequence Fidelity')
title_str = {'Sequence Fidelities from experiment',...
    [data.sets.(filenames{1}).Experiment_Name,' on ',data.sets.(filenames{1}).Timestamp]};
title(title_str)
leg = legend(h,'Data',...
    ['Mean = (',num2str(data.seq_fid_mean),'\pm',num2str(data.seq_fid_std),')%']);
set(leg,'Location','best');
end