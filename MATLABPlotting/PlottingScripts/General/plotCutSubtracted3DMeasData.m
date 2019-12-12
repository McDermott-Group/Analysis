function plotCutSubtracted3DMeasData(normalization_direction,...
    data_variable)
%plotCutSubtracted2DMeasData(NORMALIZATION_DIRECTION, DATA_VARIABLE) Plot
%a 2D data subtracting a median cut from each slice.
%   plotCutSubtracted2DMeasData(DATA_VARIABLE, NORMALIZATION_DIRECTION)
%   plots cut-substracted data. The slices are averaged and removed from
%   the data.

if exist('normalization_direction', 'var') &&...
        ~isempty(strfind(normalization_direction, 'x'))
    normalization_direction = 'along_x';
else
    normalization_direction = 'along_y';
end

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

if ~exist('data_variable', 'var')
    data_variable = selectDepDataVars(data, true);
    if isempty(data_variable)
        return
    end
    data_variable = data_variable{1};
end

% Check that the data variable exists (compute it if necessary).
[data, data_variable] = checkDataVar(data, data_variable);

[~, filename, ~] = fileparts(data.Filename);

dep_vals = data.(data_variable);
dep_rels = data.rels.(data_variable);
    
if isempty(dep_rels)
    error(['Independent (sweep) variables for data variable ''',...
          strrep(data_variable, '_', ' '), ''' are not specified.'])
end


% Plot 2D extracted frequency from 3D data
if length(dep_rels) == 3
    if strcmp(normalization_direction, 'along_y')
        
        for i = 1:size(dep_vals, 2)
            dep_vals(:,i,:) = squeeze(dep_vals(:,i,:)) - (ones(size(squeeze(dep_vals(:,i,:)),1),1) *...
                median(squeeze(dep_vals(:,i,:))));
        end 
    elseif strcmp(normalization_direction, 'along_x')
        dep_vals = dep_vals - (median(dep_vals, 2) *...
            ones(1, size(dep_vals, 2)));
    end
        
    for i = 1:size(dep_vals, 1)
        for j = 1:size(dep_vals, 2)
            [t, k] =  min(dep_vals(i,j,:));
%             if data.RF_Frequency(k) > 6.445 |  t > -0.1
%                 min_dep_vals(i,j) = 6.36;
%                 min_dep_valsS21(i,j) = -0.05;
%             else 
                min_dep_vals(i,j) = data.RF_Frequency(k)
                min_dep_valsS21(i,j) = t;
%             end
        end 
    end
    
    processed_data_var = ['minimumOfCutSubtracted_', dep_rels{3}];
    processed_data_var2 = ['minimumOfCutSubtracted_', data_variable];
    data.(processed_data_var) = min_dep_vals;
    data.(processed_data_var2) = min_dep_valsS21;
    data.units.(processed_data_var) = data.units.RF_Frequency;
    data.rels.(processed_data_var) = data.rels.(data_variable);
    data.units.(processed_data_var2) = data.units.(data_variable);
    data.rels.(processed_data_var2) = data.rels.(data_variable);
    data.dep{length(data.dep) + 1} = processed_data_var;
    data.dep{length(data.dep) + 1} = processed_data_var2;
    data.plotting.(processed_data_var).full_name =...
        ['min-Cut-Subtracted ', strrep( dep_rels{3}, '_', ' ')];
    data.plotting.(processed_data_var).extra_filename =...
        ['_', normalization_direction];

    %Uncomment if you want to plot extracted minimum from 3D dataset
    plotDataVar(data, processed_data_var);
    %colormap(winter)
  
    
    
    %For looking at cuts from 3D dataset
%     figure;
%     h = plot(data.(dep_rels{3}), squeeze(dep_vals(:,1,:)), '-', 'LineWidth',2);
%     xlabel('RF Frequency(GHz)');
%     ylabel('S21(dB)');
%     title('Cuts from 3D dataset', 'FontSize', 10);
%     grid on;
%     axis tight;
  
    
    saveMeasData(data, [filename, '_', data_variable, 'min_cut_subtr'])
end

