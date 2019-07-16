classdef chargeScan
    properties
        qubit
        file_tag
        start = 1;
        stop = 0;
        measurement_time
        psd
        f
        time
        raw_voltage
        unwrapped_voltage
        dedv
        wrapping_voltage
        unwrap_mult
        % properties for each of the four species of qubits:
        dedv_list = [0.3844, 0.411, nan, 0.5754];
        wrapping_voltage_list = [0.3844, 0.411, nan, 0.5754];
    end
    methods
        function o = chargeScan(file_tag, qubit, start, stop)
            o.file_tag = file_tag;
            o.qubit = qubit;
            o.start = start;
            o.stop = stop;
            fileName = strcat(['Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-06-10\Circ1\Q', ...
                                num2str(qubit), '\General\Parameter\', ...
                                file_tag, '_parameters.hdf5']);
            charge_file = loadMeasurementData(fileName);
            
            if (stop==0)
                o.time = charge_file.Time(start:end);
                o.raw_voltage = charge_file.Offset_Voltage(start:end) - charge_file.Offset_Voltage(1);
            else
                o.time = charge_file.Time(start:stop);
                o.raw_voltage = charge_file.Offset_Voltage(start:stop) - charge_file.Offset_Voltage(1);
            end
            
            o.dedv = 2/o.dedv_list(qubit);%e/V
            o.wrapping_voltage = o.wrapping_voltage_list(qubit)/2;%V
            o.unwrap_mult = 1;%If we unwrap at 5 V, then we are dropping 2e.  dedv converts
            %voltage to 1e, so we need a multiplier for the unwrapping.
            
            o.measurement_time = (o.time(end) - o.time(1))/length(o.time);
            l = length(o.raw_voltage);
            voltages = o.raw_voltage(1:l - mod(l,2));
            [o.psd,o.f,o.unwrapped_voltage,t1] = Wrapper_AnalyzeChargeTracePSD(voltages,o.wrapping_voltage,o.unwrap_mult,...
                            o.dedv,o.measurement_time);
        end
        function o = find_breaks(o)
            T = o.time(2:end) - o.time(1:end-1);
            breaks = find(T>5*o.measurement_time)
            nBreaksForFill = (o.time(breaks+1)-o.time(breaks))/o.measurement_time
        end
        function o = scan_length(o)
            scan_length = (o.time(end)-o.time(1))/60/60
        end
        function o = plot_psd(o)
            figure(150); hold on;
            title('PSDs')
            plot(o.f,o.psd, 'DisplayName', strcat(['Q',num2str(o.qubit)]))
            set(gca,'xscale','log')
            set(gca,'yscale','log')
            xlabel('Frequency (Hz)')
            ylabel('S_q (e^2/Hz)')
            grid on
        end
        function [f, psd] = calc_psd(o, l, measurement_time)
            [psd,f,unwrapped_voltage,t1] = Wrapper_AnalyzeChargeTracePSD(o.raw_voltage(1:l),o.wrapping_voltage,o.unwrap_mult,...
                            o.dedv,measurement_time);
        end
        function o = show_histogram(o)
            [delta_1] = calc_delta(o.unwrapped_voltage,1);
            delta_1 = mod(0.5 + delta_1, 1) - 0.5;
            figure(130+o.qubit); hold on;
            h=histogram(delta_1,50);
            edges = h.BinEdges;
            counts = h.BinCounts;
            totalTime = o.time(end) - o.time(1);
            histogram('BinEdges', edges, 'BinCounts', counts/totalTime)
            xlabel('Jump size (e)')
            ylabel('N')
            set(gca,'yscale','log')
        end
    end
end



function [delta] = calc_delta(es,stepsize)
    delta = zeros(length(es)-stepsize,1);
    for i = 1:length(delta)
        delta(i) = es(i+stepsize) - es(i);
    end
end

function [es_jump, es_smooth] = filterJumps(es,delta)
thresh = 0.07;%0.001;
es_jump = zeros(length(es),1) + es(1);
es_smooth = zeros(length(es),1) + es(1);
for i = 1:length(delta)
    if abs(delta(i)) < thresh
        es_smooth(i+1) = es_smooth(i) + delta(i);
        es_jump(i+1) = es_jump(i);
    else
        es_jump(i+1) = es_jump(i) + delta(i);
        es_smooth(i+1) = es_smooth(i);
    end
end
end

function [psd_cut,freqs] = crosspsd(z1,z2,ts)
Fs = length(ts)/(ts(end) - ts(1));
NumP = length(z1);
freqs = 0:Fs/NumP:(Fs/2);
%Cross PSD
NumP = length(z1);
fft_seq1 = fft(z1);
fft_seq2 = conj(fft(z2));
psd = (1/(Fs*NumP))*(fft_seq1.*fft_seq2);
psd(2:(NumP/2)) = 2*psd(2:(NumP/2));
%End PSD
psd_cut = psd(1:(NumP/2+1));
end

function [apsd] = window_averaging(psd)
% apsd = psd;
%Averaging with f window
apsd = zeros(size(psd));
for i = 1:size(psd,1)
    filter_fl = max(1,round(i - i/4));
    filter_fh = min(size(psd,1),round(i + i/4));
    apsd(i) = mean(psd(filter_fl:filter_fh));
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