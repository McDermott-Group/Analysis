function fitMeasData2avoidedcrossing_ameya(data_variable)

% Select a file.
if ~exist('data', 'var')
    % Select a file.
    data = loadMeasurementData;
end
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

[pathname, filename, ext] = fileparts(data.Filename);

% Create folder Plots if necessary.
plts_path = makeDirPlots(pathname);

% Check that the data variable exists (compute it if necessary).
[data, data_variable] = checkDataVar(data, data_variable);

dep_vals = data.(data_variable);
dep_rels = data.rels.(data_variable);

if isempty(dep_rels)
    error(['Independent (sweep) variables for data variable ''',...
          strrep(data_variable, '_', ' '), ''' are not specified.'])
end

if ~isempty(strfind(dep_rels{1}, 'JPM_Flux_2')) ||...
        ~isempty(strfind(dep_rels{1}, 'RF_Frequency')) ||...
        ~isempty(strfind(dep_rels{1}, 'Qubit_Drive_to_Readout'))
    indep_name1 = dep_rels{2};
    indep_name2 = dep_rels{1};
    indep_vals1 = data.(dep_rels{2});
    indep_vals2 = data.(dep_rels{1});
    dep_vals = dep_vals';
    flip_fit = true;
elseif ~isempty(strfind(dep_rels{2}, 'JPM_Flux_2')) ||...
        ~isempty(strfind(dep_rels{2}, 'RF_Frequency')) ||...
        ~isempty(strfind(dep_rels{2}, 'Qubit_Drive_to_Readout'))
    indep_name1 = dep_rels{1};
    indep_name2 = dep_rels{2};
    indep_vals1 = data.(dep_rels{1});
    indep_vals2 = data.(dep_rels{2});
    flip_fit = false;
else
    error(['The data does not appear to depenend on any ',...
        '''Time'', ''Duration'' or ''Qubit Drive to Readout'''...
        'variables.'])
end
dep_units = data.units.(data_variable);
indep2_units = data.units.(indep_name2);

%plot 2d data
if ~exist('data_variable', 'var')
    dep_vars = selectDepDataVars(data);
    if isempty(dep_vars)
        return
    end
    for data_index = 1:length(dep_vars)
        dep_name = dep_vars{data_index};
        if contains(dep_name, '_Std_Dev') ||...
                 contains(dep_name, '_Error')
            continue
        end
        if contains(dep_name, 'Phase')
            data.(dep_name) = data.(dep_name); %unwrap(data.(dep_name));
        end
        plotDataVar(data, dep_name);
    end
else
    plotDataVar(data, data_variable);
    hold on
end
temp_g=1e9;
minima=zeros(1,2);
fluxbias=zeros(1,1);
for j=1:length(indep_vals2)
    TF1 = islocalmin(dep_vals(:,j),'MinProminence',2);
    TF=find(TF1);
%     plot(indep_vals1, dep_vals(:,15), indep_vals1(TF), dep_vals(TF,15),'r*')
%     fprintf('%f',length(TF))
%     pause
    if length(TF)==2
%         fprintf('2 minima found! \n')
        minima(j,:)=indep_vals1(TF);
        fluxbias(j)=indep_vals2(j);
        if abs(minima(j,1)-minima(j,2))<temp_g
            temp_g=abs(minima(j,1)-minima(j,2));
        end
    end
end

scatter(fluxbias,minima(:,1),'filled')
hold on
scatter(fluxbias,minima(:,2),'filled')
hold off
fprintf('2g=%f \n',temp_g)

