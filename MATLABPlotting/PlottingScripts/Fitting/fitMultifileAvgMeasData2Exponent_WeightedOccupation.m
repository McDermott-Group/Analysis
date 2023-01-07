function fitMultifileAvgMeasData2Exponent()
% Select multiple 1-D Exponential sweeps for input. This script averages them and fit
% them to find decay rate. Error bars will be shown 


% Select files to plot.
[filenames, pathnames, status] = selectMeasurementDataFile;
if ~status
    return
end

if ~iscell(filenames)
    filenames = {filenames};
end
if ~iscell(pathnames)
    pathnames = {pathnames};
end


% Read the data file, convert the variable names, and specify the units.
try
    data = cell(0);
    for k = 1:length(filenames)
        file = fullfile(pathnames{k}, filenames{k});
        data{k} = processMeasurementData(importMeasurementData(file));
    end
catch
    error('The selected files contain unmatched data.')
end


% if exist('errorss','var')
% errors = errorss;
% end


for k = 1:length(filenames)
    Occupation(k,:) = data{k}.Weighted_Occupation;
end


AvgData.Occupation = mean(Occupation,1);


 
dep_rels = data{1}.rels.Weighted_Occupation;

% Analyze and plot 1D data
if length(dep_rels) == 1

    AverageOccupation = AvgData.Occupation;
    Time = data{1}.QB_Idle_Gate_Time';
    
    % Fit the two curves for depolarizing parameter
    
    [out.Time, out.AverageOccupation]=prepareCurveData(Time, AverageOccupation);
    [out.frExpFit, ~] = fidelityFit(out.Time, out.AverageOccupation);
    
    
    conf = confint(out.frExpFit, 0.6827);
    
    std_a = (conf(2,1) - conf(1,1))/2;
    std_1byb = max([abs(1 /conf(2,2) - 1 / out.frExpFit.b), abs(1 / conf(1,2) - 1 / out.frExpFit.b)]);
    std_c = (conf(2,3) - conf(1,3))/2;
    

    


% r_c_est_forjustInt = 0.5*(1-out.fr_Interleave.b/b_for_Interleave_forjustInt);
% fidelity_gate_forjustInt = 1 - r_c_est_forjustInt;
% term1_forJustInt = std_i_b/b_for_Interleave_forjustInt;
% term2_forJustInt = out.fr_Interleave.b/b_for_Interleave_forjustInt^2 * std_no_b_forjustInt;
% sigma_r_c_est_justforInt = 0.5 * sqrt(term1_forJustInt^2 + term2_forJustInt^2);
% 

    
    % Plotting
    createFigure; hold on;
    
   
    %plot the fits on top of the data
    
    ax=gca;
    ax.ColorOrderIndex=1;
    h = plot(out.frExpFit, out.Time, out.AverageOccupation);
    h(1).LineStyle = 'none'; 
    h(1).Marker = 'o'; 
    h(1).MarkerSize = 6;
    h(1).MarkerFaceColor = [0, 0.447 , 0.741];
    h(1).MarkerEdgeColor = 'none';
    h(2).LineWidth = 2;
    h(2).Color = 'r';
    
    [~,filename,ext] = fileparts(data{1}.Filename);
    filename=strrep(filename,'_','\_');
    freepars = ['a = ', num2str(1/out.frExpFit.a),'±',num2str(std_a),';',...
                '1/b = ', num2str(1/out.frExpFit.b),'±',num2str(std_1byb),'nsec ;',...
                'c = ', num2str(out.frExpFit.c),'±', num2str(std_c)];
    eqn = strcat([strrep('Occupation', '_', ' '),' = a * exp(-b * ',...
                strrep('Idle Gate Time', '_', ' '), ') + c'],';');
    title({[strrep(filenames{1}, '_', '\_'), ' - ',...
    strrep(filenames{end}, '_', '\_'),': ',data{1}.Timestamp],...
        eqn, freepars}, 'FontSize', 12); 
end 

    xlabel('Idle Gate Time (nsec)')
    ylabel('Average Occupation')
    legend('data', 'fit')
    axis tight
    grid on
%     leg=legend('Interleaved Gate: None',...
%         ['Interleaved Gate: ',data{1}.Interleaving_Gate]);
%     set(leg,'Location','best')
    
end


function [fitresult, gof] = fidelityFit(x, y)

% Set up fittype and options.
    opts = fitoptions('Method', 'NonlinearLeastSquares',...
                      'Robust', 'LAR',...
                      'MaxFunEvals',10000,...
                      'MaxIter',10000,...
                      'DiffMaxChange',1,...
                      'Algorithm','Trust-Region');
    opts.Display = 'Off';
    opts.TolX = 1e-9;
    opts.TolFun = 1e-9;
    opts.StartPoint = [y(1) - y(end), 1 / (max(x) - min(x) + 9 * eps),y(end)];
    [fitresult, gof] = fit(x(:), y(:), 'a * exp(-b * x) + c', opts);

end