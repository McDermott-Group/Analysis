function [window_avg_psd, psd_freq] = Plot_QP_Tunneling(CDdate, samples, ...
                                    qubit, date, fileIndicies)

    %PARAMETERS
%     reps = 8192;
%     trials = 5;
    n_avg = 1; % number of samples to average together before taking PSD
    refreshTime = n_avg*100e-6;

    dataType = 'QP_Tunneling_PSD';
    path = strcat(['Z:\mcdermott-group\data\fluxNoise\',CDdate,samples, ...
              qubit,'General\',date,dataType,'\MATLABData']);

    Fs = 1/refreshTime; %samples per second
    vis = 1;
    transfer = 1/vis^2;%1.485e-5/0.0002957;
    
    % load the data into the needed matrix
    reps = 1; trials = 1; data = 0;
    j = 0;
    for i = fileIndicies
        ldata = noiselib.load_file(path, [dataType '_' num2str(i,'%03d') '.mat']);
        o = ldata.Single_Shot_Occupations;
        oo = zeros(5,8192/n_avg);
        if n_avg > 1
            for k = 1:5
                oo(k,:) = mean(reshape(o(k,:), n_avg, 8192/n_avg)) > 0.1;
            end
            o = oo;
%         o = movmean(o,5,2)>0.1;
        end
        if j == 0
            trials = size(o,1);
            reps = size(o,2);
            data = zeros(trials*length(fileIndicies), reps);
        end
        data(j*trials + (1:trials), :) = o;
        j = j + 1;
    end
    
    [avg_cpsd, psd_freq] = noiselib.partition_and_avg_psd(data, Fs);
    avg_cpsd = avg_cpsd' * transfer;
    window_avg_psd = noiselib.window_averaging(avg_cpsd);

    %Fit to find QP tunneling rate:
    % x = psd_freq(2:end);
    % y = abs(window_avg_psd(2:end));
    % [x, y] = prepareCurveData(x, y);
    % [fr,~] = noiselib.fit_lorenzian(x,y);

    % f = 0.01:0.01:5e3; gamma = 1056;lor = 3.8e2./(4*pi^2*f.^2+gamma^2);
    % figure(12);hold on;plot(f,lor)

    %Plot
    figure(113);hold on
    title('1/f Averaged PSD')
    plot(psd_freq,abs(window_avg_psd), 'DisplayName', [samples(1:end-1),qubit(1:end-1)])
    % plot(fr, x, y);
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    xlabel('Frequency (Hz)')
    ylabel('S_\eta (\eta^2/Hz)')
    grid on

end