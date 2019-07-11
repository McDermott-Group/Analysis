

scans = {...
    % chargeScan('ckv0334jpp', 2, 1, 865), ... % longer measurement time
    % chargeScan('ckv0334jpp', 2, 866, 0), ...
    % chargeScan('cla0423hwf', 1, 1, 1770), ... % bad fits throughout
    % chargeScan('cla0423hwf', 1, 1771, 0), ...
    % chargeScan('clf2154lqz', 1, 1, 0), ... % maybe bad fits?
    chargeScan('ckx2109xfx', 4, 182, 2402), ...
    chargeScan('ckx2109xfx', 4, 2403, 0), ...
    chargeScan('cky2155mll', 1, 1, 1674), ... % cut off before fill before bad data
    chargeScan('clb0821qqz', 4, 1, 1048), ...
    chargeScan('clb0821qqz', 4, 1048, 0), ...
    chargeScan('clc0632vzv', 2, 1, 5033), ... % cut off end when power outage happened
    chargeScan('cle0301xnh', 2, 1, 0), ...
    chargeScan('clg1638iiy', 2, 1, 0), ...
    chargeScan('clh2255qso', 4, 1, 0), ...
    chargeScan('cli1847ppc', 4, 1, 0), ...
    chargeScan('clm2340vve', 2, 1, 0), ...
    };

% scans{6}.find_breaks()

for q = [1 2 4]
    % combine voltages and plot all unwrapped jumps for given qubit
    combined_voltages = [0];
    combined_time = [0];
    meas_times = [];
    scan_lengths = [];
    for i = 1:length(scans)
        if scans{i}.qubit == q
            combined_voltages = [combined_voltages, combined_voltages(end) + scans{i}.unwrapped_voltage'];
            l = length(scans{i}.unwrapped_voltage);
            combined_time = [combined_time, scans{i}.time(1:l)'-scans{i}.time(1)+combined_time(end)];
            meas_times = [meas_times, scans{i}.measurement_time];
            scan_lengths = [scan_lengths, length(scans{i}.unwrapped_voltage)];
        end
    end
    figure(200); hold on;
    plot(combined_time, combined_voltages, 'DisplayName', strcat(['Q',num2str(q)]))
    legend;
    
    % plot histogram and histogram normalized by total time
    figure(210); hold on;
    [delta_1] = calc_delta(combined_voltages,1);
    delta_1 = mod(0.5 + delta_1, 1) - 0.5;
    h=histogram(delta_1,50, 'DisplayName', strcat(['Q',num2str(q)]));
    xlabel('Jump size (e)')
    ylabel('N')
    set(gca,'yscale','log')
    legend;
    edges = h.BinEdges;
    counts = h.BinCounts;
    totalTime = combined_time(end) - combined_time(1);
    figure(211); hold on;
    histogram('BinEdges', edges, 'BinCounts', counts/totalTime, 'DisplayName', strcat(['Q',num2str(q)]));
    xlabel('Jump size (e)')
    ylabel('N')
    set(gca,'yscale','log')
    legend;
    sum(abs(delta_1) > 0.1)/totalTime
    
    % plot averaged psds
    avg_psd = 0;
    for i = 1:length(scans)
        if scans{i}.qubit == q
            scans{i}.plot_psd();
            l = min(scan_lengths) - mod(min(scan_lengths), 2); %lowest even number
            [f, psd] = scans{i}.calc_psd(l, mean(meas_times));
            figure(150); hold on;
            plot(f, psd)
            if avg_psd == 0
                avg_psd = psd/length(meas_times);
            else
                avg_psd = avg_psd + psd/length(meas_times);
            end
        end
    end
    figure(150); hold on;
    title('PSDs')
    plot(f,avg_psd, 'DisplayName', strcat(['Q',num2str(q), ' Averaged']))
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel('S_q (e^2/Hz)')
    grid on
    
    % fit
    [fr,~] = fit_psd(log(f),log(avg_psd));
    alpha = fr.b;
    A = fr.a;
    [alpha, A]
    fVals = [f(1), f(end), 1e0];
    plot(fVals, exp(A + log(fVals)*alpha), 'DisplayName', strcat(['Q',num2str(q), ' Fit']))
end

function [delta] = calc_delta(es,stepsize)
    delta = zeros(length(es)-stepsize,1);
    for i = 1:length(delta)
        delta(i) = es(i+stepsize) - es(i);
    end
end

function [fitresult, gof] = fit_psd(f,psd)
[xData, yData] = prepareCurveData(f, psd');

ft = fittype( 'a + b*x', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [-7 -1.5];

[fitresult, gof] = fit( xData, yData, ft, opts );

% figure;plot( fitresult, xData, yData );
end
