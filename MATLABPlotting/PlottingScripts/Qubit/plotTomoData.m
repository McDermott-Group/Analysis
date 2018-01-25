function [pax,cbar]=plotTomoData(loadedTomoData)
%plotTomoData   Plot tomography data from a data file.
%   plotTomoData plots tomography data from a selected data file.

% Additional functionality Oct 2016: able to specify an already loaded data
%   structure that has been loaded via loadMeasurementData

% Select a file if loadedTomoData doesn't exist
if ~exist('loadedTomoData', 'var')
    data = loadMeasurementData;
else
    data = loadedTomoData;
end

if isempty(fields(data))
    return
end

data_variables = {'Sigma_X','Sigma_Y','Sigma_Z'};

% Check that the tomography data variables exist and are all 1D arrays.
for k = 1:length(data_variables)
    if ~isfield(data, data_variables{k})
        error('Data file does not contain complete tomography data.')
    end
    if length(data.rels.(data_variables{k})) ~= 1
        error(['Data for variable ''', data_variables{k}, ''' is not 1D.'])
    end
end

[pathname, filename, ext] = fileparts(data.Filename);
plts_path = makeDirPlots(pathname);
plot_title = [filename, ext, ' [', data.Timestamp, ']'];

indep = data.rels.(data_variables{1}){1};
xunits = getUnits(data, indep);
xlables = [strrep(indep, '_', ' '), xunits];
yunits = getUnits(data, data_variables{1});
ylables = data_variables;
legends = data_variables;

for k = 1:length(data_variables)
    if ~strcmp(indep, data.rels.(data_variables{k}){1})
        error(['The independent variables are not the same ',...
            'for specified data variables.'])
    end
    units = getUnits(data, data_variables{k});
    if ~strcmp(yunits, units)
        error('The dependent variable units do not match.')
    end
    ylables{k} = [strrep(data_variables{k}, '_', ' '), units];
    legends{k} = strrep(data_variables{k}, '_', ' ');
end

% Plot the expectation values
createFigure;
hold on
for k = 1:length(data_variables)
    if strcmp(data_variables{k},'Phase')
        data.(data_variables{k}) = unwrap(data.(data_variables{k}));
    end
   
    indeps = data.rels.(data_variables{k});
    plot(data.(indeps{1}), data.(data_variables{k}),...
        '.-', 'LineWidth', 1, 'MarkerSize', 15)
end
hold off
axis tight
grid on
set(gca, 'box', 'on')

xlabel(xlables, 'FontSize', 14)
ylabel(ylables, 'FontSize', 14)
title(plot_title, 'Interpreter', 'none', 'FontSize', 10)
legend(legends, 'Interpreter', 'none', 'Location', 'Best')
savePlot(fullfile(plts_path, [filename, '_comb_simple']));

% plot the XZ plane
figure;
subplot(1,2,1)
th_xz = atan2(-data.Sigma_Z,data.Sigma_X);
r_xz = sqrt(data.Sigma_Z.^2 + data.Sigma_X.^2);
pax=polarscatter(th_xz,r_xz,[],data.(indeps{1}),'filled');
colormap jet
cbar=colorbar;
cbar.Label.String = xlables;
cbar.Label.Rotation = 270;
cbar.FontSize=14;
cbar.Label.Position=cbar.Label.Position + [1 0 0];
thetaticks([0 90 180 270])
thetaticklabels({'+X','|0> (-Z)','-X','|1> (+Z)'})
title(plot_title, 'Interpreter', 'none', 'FontSize', 10)

% plot the YZ plane

subplot(1,2,2)
th_yz = atan2(-data.Sigma_Z,data.Sigma_Y);
r_yz = sqrt(data.Sigma_Z.^2 + data.Sigma_Y.^2);
pax=polarscatter(th_yz,r_yz,[],data.(indeps{1}),'filled');
colormap jet
cbar=colorbar;
cbar.Label.String = xlables;
cbar.Label.Rotation = 270;
cbar.FontSize=14;
cbar.Label.Position=cbar.Label.Position + [1 0 0];
thetaticks([0 90 180 270])
thetaticklabels({'+Y','|0> (-Z)','-Y','|1> (+Z)'})
title(plot_title, 'Interpreter', 'none', 'FontSize', 10)
savePlot(fullfile(plts_path, [filename, '_tomo_YZ']));


