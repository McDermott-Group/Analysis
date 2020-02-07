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
        raw_time
        raw_voltage
        unwrapped_voltage
        dedv
        wrapping_voltage
    end
    methods
        function o = chargeScan(file_tag, qubit, start, stop)
            o.file_tag = file_tag;
            o.qubit = qubit;
            o.start = start;
            o.stop = stop;
            fileName = strcat(['Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-08-12\Circ1\Q', ...
                                num2str(qubit), '\General\Parameter\', ...
                                file_tag, '_parameters.hdf5']);
            charge_file = loadMeasurementData(fileName);
            
            offset = 0.;
            
            if (stop==0)
                o.raw_time = charge_file.Time(start:end);
                o.raw_voltage = charge_file.Offset_Voltage(start:end)-offset;
            else
                o.raw_time = charge_file.Time(start:stop);
                o.raw_voltage = charge_file.Offset_Voltage(start:stop)-offset;
            end
            
            o.dedv = 2/charge_file.x2e_Period; %e/V
            o.wrapping_voltage = charge_file.x2e_Period/2;%V
            
            o.measurement_time = (o.raw_time(end) - o.raw_time(1))/length(o.raw_time);
            l = length(o.raw_voltage);
            voltages = o.raw_voltage(1:l - mod(l,2));
            [o.psd,o.f,o.unwrapped_voltage,o.time] = o.analyze_voltage_trace(voltages,o.wrapping_voltage,...
                            o.dedv,o.measurement_time);
            o.unwrapped_voltage = o.unwrapped_voltage - o.unwrapped_voltage(1);
        end
        function o = find_breaks(o)
            T = o.raw_time(2:end) - o.raw_time(1:end-1);
            breaks = find(T>5*o.measurement_time)
            nBreaksForFill = (o.raw_time(breaks+1)-o.raw_time(breaks))/o.measurement_time
        end
        function o = scan_length(o)
            scan_length = (o.raw_time(end)-o.raw_time(1))/60/60
        end
        function o = plot_unwrapped_voltages(o)
            figure(140); hold on;
            title('Voltages')
            plot(o.time - o.time(1),o.unwrapped_voltage, 'DisplayName', strcat(['Q',num2str(o.qubit)]))
            xlabel('Time (s)')
            ylabel('Voltage (e)')
            grid on
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
            [psd,f,unwrapped_voltage,t1] = o.analyze_voltage_trace(o.raw_voltage(1:l),o.wrapping_voltage,...
                            o.dedv,measurement_time);
        end
        function [aqpsd,f,es,t] = analyze_voltage_trace(vs,wrapping_voltage,...
                                         dedv, measurement_time)
            Fs = 1/measurement_time;%samples per second
            NumP = length(vs);
            f = 0:Fs/NumP:(Fs/2);

            [es] = noiselib.unwrap_voltage_to_charge(vs, wrapping_voltage, dedv);
            t = (0:(length(vs)-1))*(measurement_time);

            [Qpsd] = noiselib.crosspsd(es,es,Fs/2);
            [aqpsd] = noiselib.window_averaging(Qpsd);
            aqpsd = aqpsd(1:(NumP/2+1));
        end
        function o = show_histogram(o)
            [delta_1] = noiselib.calc_delta(o.unwrapped_voltage,1);
            delta_1 = noiselib.alias(delta_1, 0.5);
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