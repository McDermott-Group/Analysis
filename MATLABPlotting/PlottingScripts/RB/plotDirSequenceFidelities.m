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
        fprintf(fid, '%s', dataDir);
        fclose(fid);
    end
else
    data = loadedDirData;
end

% With the data loaded, grab the occupations and gate sequences of each run
filenames = fieldnames(data.sets);
data.seq_fids = zeros(length(filenames),3);
data.seq_times = zeros(length(filenames),1);
data.seq_lengths = zeros(length(filenames),1);
data.gate_seqs{length(filenames),1} = [];


data.Pi_Time = data.sets.(filenames{1}).SFQ_RB_Pulse_Duration;
data.Pi2_Time = data.sets.(filenames{1}).SFQ_RB_Pi_Over_2_Duration;

data.gates = {'I','X','Y','X2','Y2'};
data.gate_lengths = [0,1,1,0.5,0.5];
data.gate_times = [0, data.Pi_Time, data.Pi_Time, data.Pi2_Time,...
                      data.Pi2_Time];

for k=1:length(filenames)
   data.seq_fids(k,1) = 1 - data.sets.(filenames{k}).Poor_Mans_Occupation;
   data.gate_seqs{k} = regexprep(...
       data.sets.(filenames{k}).SFQ_User_Defined_Sequence,']|[|''| ','');
   simple_seq = regexprep(data.gate_seqs{k},'-|/','');
   gates=strsplit(simple_seq,',');
   num_gates = length(gates);
   for gate_idx = 1:num_gates
      for check_idx = 1:length(data.gates)
         if gates{gate_idx} == data.gates{check_idx}
            data.seq_times(k) = data.seq_times(k) + data.gate_times(check_idx);
            data.seq_lengths(k) = data.seq_lengths(k) + data.gate_lengths(check_idx);
         end
      end
   end
   if data.seq_lengths(k) == 0
      data.measure_fidelity = data.seq_fids(k,1); 
   end
end

% normalize one set of data by the measurement fidelity
data.seq_fids(:,2) = data.seq_fids(:,1)./data.measure_fidelity;

% % normalize another set by sequence length
% for j = 1:4
%    data.gate_length_idx(j,:) = data.seq_lengths==j;
%    data.gate_length_norm_fids(j) = mean(data.seq_fids(data.gate_length_idx(j,:),1));
%    data.gate_adjustment_factor(j) = 1;
%    if j>1
%       data.gate_adjustment_factor(j) = data.gate_length_norm_fids(1) / ...
%                                        data.gate_length_norm_fids(j);
%    end
% end
% 
% for j=1:length(data.seq_fids)
%     if data.seq_lengths(j) == 0
%         data.seq_fids(j,3) = data.seq_fids(j,2);
%     else
%         data.seq_fids(j,3) = data.seq_fids(j,2) * data.gate_adjustment_factor(data.seq_lengths(j));
%     end
% end

% divide out the contributions from the Rabi envelope (hardcoded for now)
data.rabi_envelope_time = 1673; % ns
for j=1:length(data.seq_fids)
    data.seq_fids(j,3) = data.seq_fids(j,1) / exp(-data.seq_times(j)/(2*data.rabi_envelope_time));
end

% now plot it with its mean and std on the same plot
createFigure; hold on;
h(1:2)=bar([data.seq_fids(:,1), data.seq_fids(:,2)],'FaceColor','flat');
data.seq_fid_mean = mean(data.seq_fids(:,1));
data.adj_seq_fid_mean = mean(data.seq_fids(:,2));
data.seq_fid_std = std(data.seq_fids(:,1));
data.adj_seq_fid_std = std(data.seq_fids(:,2));
h(3)=refline(0,data.seq_fid_mean);
h(3).Color = 'k';
h(4)=refline(0,data.adj_seq_fid_mean);
h(4).Color = 'b';
% hline=refline(0,data.seq_fid_mean+data.seq_fid_std);
% hline.Color = 'r';
% hline=refline(0,data.seq_fid_mean-data.seq_fid_std);
% hline.Color = 'r';
ax=gca;
set(ax,'XTickLabelRotation',300)
set(ax,'XTick',1:1:length(filenames))
set(ax,'XTickLabel',data.gate_seqs)
set(ax,'YGrid','on')
set(ax,'YMinorGrid','on')
set(ax,'XGrid','off')
set(ax,'FontSize',12)
%set(ax,'YLim',[0,1])
xlabel('Gate Sequence')
ylabel('Sequence Fidelity')
title_str = {'Sequence Fidelities from experiment',...
    [data.sets.(filenames{1}).Experiment_Name,' on ',data.sets.(filenames{1}).Timestamp]};
title(title_str)
leg = legend(h,'Raw Data',...
    'Adjusted for Identity Infidelity',...
    ['Raw Mean = (',num2str(data.seq_fid_mean),'\pm',num2str(data.seq_fid_std),')%'],...
    ['Adjusted Mean = (',num2str(data.adj_seq_fid_mean),'\pm',num2str(data.adj_seq_fid_std),')%']);
set(leg,'Location','best');

createFigure; hold on;
nbins = 13;
histfit(data.seq_fids(:,2),nbins)
xlabel('Adjusted Sequence Fidelity')
ylabel('Counts')
data.hist_info = fitdist(data.seq_fids(:,2),'Normal');

end