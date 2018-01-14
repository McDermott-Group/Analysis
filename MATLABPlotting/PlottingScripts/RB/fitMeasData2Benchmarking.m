function [out]=fitMeasData2Benchmarking(loadedRBData)

% Select a file if loadedRBData doesn't exist
if ~exist('loadedRBData', 'var')
    data = loadMeasurementData;
else
    data = loadedRBData;
end

% Check loaded data for possibility of generating errorbars
if isfield(data, 'Interleaved_Probabilities')
    errors = 1;
else
    errors = 0;
end

dep_rels = data.rels.Interleaved_Probability;

% analyze and plot 1D data
if length(dep_rels) == 1
    %Find P1 for Interleaved and no interleave for each number of Cliffords
    p1_array = zeros(length(data.Number_of_Cliffords), 3);%NoInter, Inter, NoC
    p1_array(:,1) = 1 - data.No_Interleaved_Probability;
    p1_array(:,2) = 1 - data.Interleaved_Probability;
    p1_array(:,3) = data.Number_of_Cliffords;
    
    % Fit the two curves for depolarizing parameter
    
    [out.no_int_x, out.no_int_y]=prepareCurveData(p1_array(:,3), p1_array(:,1));
    [out.fr_No_Interleave, ~] = fidelityFit(out.no_int_x, out.no_int_y);
    
    [out.int_x, out.int_y]=prepareCurveData(p1_array(:,3), p1_array(:,2));
    [out.fr_Interleave, ~] = fidelityFit(out.int_x, out.int_y);
    
    r_c_est = 0.5*(1-out.fr_Interleave.b/out.fr_No_Interleave.b);
    out.fidelity_gate = 1 - r_c_est;
    
    fr_i_conf = confint(out.fr_Interleave);
    fr_no_conf = confint(out.fr_No_Interleave);
    std_i_b = (fr_i_conf(2,2) - fr_i_conf(1,2))/4;
    std_no_b = (fr_no_conf(2,2) - fr_no_conf(1,2))/4;
    
    term1 = std_i_b/out.fr_No_Interleave.b;
    term2 = out.fr_Interleave.b/out.fr_No_Interleave.b^2 * std_no_b;
    out.sigma_r_c_est = 0.5 * sqrt(term1^2 + term2^2);
    
    % Plotting
    createFigure; hold on;
    
    %plot data with or without errorbars based on whether or not the datafile
    % contains the required information
    if errors
        out.no_int_std = std(data.No_Interleaved_Probabilities,1,2);
        out.int_std = std(data.Interleaved_Probabilities,1,2);
        plotErrorbar(out.no_int_x,out.no_int_y,out.no_int_std);
        plotErrorbar(out.int_x,out.int_y,out.int_std);
    else
        plotSimple(out.no_int_x,out.no_int_y);
        plotSimple(out.int_x,out.int_y);
    end
    
    %plot the fits on top of the data
    
    ax=gca;
    ax.ColorOrderIndex=1;
    plotSimple(data.Number_of_Cliffords,...
        out.fr_No_Interleave(data.Number_of_Cliffords),'-');
    plotSimple(data.Number_of_Cliffords,...
        out.fr_Interleave(data.Number_of_Cliffords),'-');
    
    [~,filename,ext] = fileparts(data.Filename);
    filename=strrep(filename,'_','\_');
    fid_str = [data.Interleaving_Gate,' Gate Fidelity: ',...
        num2str(out.fidelity_gate),'\pm',num2str(out.sigma_r_c_est)];
    title({[filename,ext,': ',data.Timestamp],...
        fid_str});
    xlabel('Number of Cliffords')
    ylabel('Sequence Fidelity')
    leg=legend('Interleaved Gate: None',...
        ['Interleaved Gate: ',data.Interleaving_Gate]);
    set(leg,'Location','best')
    
% analyze and plot 2D data
elseif length(dep_rels) == 2
    
    % figure out what the swept variable is aside from number of cliffords
    if ~isempty(strfind(dep_rels{1}, 'Cliffords'))
        indep_name1 = dep_rels{2};
        indep_name2 = dep_rels{1};
        indep_vals1 = data.(dep_rels{2});
        indep_vals2 = data.(dep_rels{1});
        dep_vals = dep_vals';
        flip_fit = true;
    % elseif ~isempty(strfind(dep_rels{1}, 'Duration')) || ...
    % elseif ~isempty(strfind(dep_rels{2}, 'Delay')) || ...
    elseif ~isempty(strfind(dep_rels{2}, 'Cliffords'))
        indep_name1 = dep_rels{1};
        indep_name2 = dep_rels{2};
        indep_vals1 = data.(dep_rels{1});
        indep_vals2 = data.(dep_rels{2});
        flip_fit = false;
    else
        error(['The data does not appear to depenend on any ',...
               '''Cliffords'' variable.'])
    end
    %Find P1 for Interleaved and no interleave
    out.p1_array(:,:,1) = 1 - data.No_Interleaved_Probability;
    out.p1_array(:,:,2) = 1 - data.Interleaved_Probability;
    out.m = indep_vals2;
    out.swept_vals = indep_vals1;
    out.swept_var = indep_name1;
    
    %% For loop to process each dataset and generate relevant plots
    out.gate_fids = zeros(length(out.swept_vals),1);
    out.fidelity_gate = zeros(length(out.swept_vals),1);
    out.sigma_r_c_est = zeros(length(out.swept_vals),1);
    out.no_int_seq_fid = zeros(length(out.swept_vals), length(out.m));
    out.int_seq_fid = zeros(length(out.swept_vals), length(out.m));
    for k=1:length(out.swept_vals)
        sweep_str = ['sweep',num2str(k)];
        
        % Fit the two curves for depolarizing parameter
        [out.m, out.no_int_seq_fid(k,:)]=prepareCurveData(out.m, squeeze(out.p1_array(k,:,1)'));
        [out.fit_No_Interleave.(sweep_str), ~] = fidelityFit(out.m, out.no_int_seq_fid(k,:)');
        
        [out.m, out.int_seq_fid(k,:)]=prepareCurveData(out.m, squeeze(out.p1_array(k,:,2)'));
        [out.fit_Interleave.(sweep_str), ~] = fidelityFit(out.m, out.int_seq_fid(k,:)');
        
        r_c_est = 0.5*(1-out.fit_Interleave.(sweep_str).b/out.fit_No_Interleave.(sweep_str).b);
        out.fidelity_gate(k) = 1 - r_c_est;

        fr_i_conf = confint(out.fit_Interleave.(sweep_str));
        fr_no_conf = confint(out.fit_No_Interleave.(sweep_str));
        std_i_b = (fr_i_conf(2,2) - fr_i_conf(1,2))/4;
        std_no_b = (fr_no_conf(2,2) - fr_no_conf(1,2))/4;

        term1 = std_i_b/out.fit_No_Interleave.(sweep_str).b;
        term2 = out.fit_Interleave.(sweep_str).b/out.fit_No_Interleave.(sweep_str).b^2 * std_no_b;
        out.sigma_r_c_est(k) = 0.5 * sqrt(term1^2 + term2^2);
        
        
    
    end
    
    % plot gate fidelities as function of swept value
    createFigure; hold on;
    plotErrorbar(out.swept_vals, out.fidelity_gate, out.sigma_r_c_est)
    
    [~,filename,ext] = fileparts(data.Filename);
    filename=strrep(filename,'_','\_');
    title({[filename,ext,': ',data.Timestamp],...
        [data.Interleaving_Gate,' Gate Fidelity']});
    xlabel([strrep(out.swept_var,'_',' '),' (', data.units.(out.swept_var),')'])
    ylabel([data.Interleaving_Gate,' Gate Fidelity'])
    
end

end

function [fitresult, gof] = fidelityFit(xData, yData)

% Set up fittype and options.
ft = fittype( 'a*b^x + c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares',...
                    'Robust', 'LAR',...
                      'MaxFunEvals',10000,...
                      'MaxIter',10000,...
                      'Algorithm','Trust-Region');
opts.Display = 'Off';
opts.Lower = [0, 0, 0];
opts.StartPoint = [0.4 0.9 mean(yData)];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data
%figure(12);hold on
%plot( fitresult, xData, yData );
end