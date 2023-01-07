% % with TWPA
% dates = {'02-28-20','02-29-20','03-02-20','03-02-20','03-03-20','03-05-20'};
% start_files = [0,  0,     0,  0*638, 0,  0];
% end_files   = [847,0*1044,637,0*1559,770,0*207];

% % no TWPA
% dates = {'03-05-20'};
% start_files = [297];
% end_files   = [1168];

% % calibrated before each measurement (?)
% dates = {'03-08-20'};
% start_files = [0];
% end_files   = [162];

% optimized readout, fit each frame to calibrate
dates = {'03-12-20','03-14-20'};
start_files = [0,0];
end_files   = [392,366];

n_files = sum(end_files - start_files) + length(start_files);
n_trials = 10;
n_reps = 10000;
f = waitbar(0, 'Progress');

X4_event_date = {};
X4_event_file = [];
X4_event_trial = [];
X4_event_rep = [];
n_4X_events = 0;
buffer_radius = 500;

qubit_list = {'q1','q2','q3','q4','q12','q34','q123','q234','q341','q412','q1234'};

autocorr_struct = struct();
hist_struct = struct();
chunk_sizes = 10*[2.5,5,10,25,50,100,200,500,1000];
bin_edges = -0.5:400.5;
for i = 1:length(qubit_list)
    qs = qubit_list{i};
    autocorr_struct.(qs) = zeros(1,2*buffer_radius+1);
    autocorr_struct.([qs '_baseline']) = zeros(1,2*buffer_radius+1);
    hist_struct.(qs) = struct();
    for r = chunk_sizes
        hist_struct.(qs).(['r' num2str(r)]) = zeros(1,length(bin_edges)-1);
    end
end

avg_1_state = zeros(4, n_files);

nf = 0;
for i = 1:length(dates)
    d = dates{i};
    pathname = ['Z:\mcdermott-group\data\fluxNoise\DR1 - 2019-12-17\CorrFar\Q1Q2Q3Q4Corr\General\' d '\QP_correlation100\MATLABData'];
    % pathname = ['/Volumes/smb/mcdermott-group/data/fluxNoise/DR1 - 2019-12-17/CorrFar/Q1Q2Corr/General/' d '/Charge_resetting/MATLABData'];
    
    for j = start_files(i):end_files(i)
        
        waitbar(nf/n_files, f, [num2str(nf) '/' num2str(n_files)]);
        
        nf = nf + 1;
        filename = ['QP_correlation100_',num2str(j,'%03d'),'.mat'];
        data = noiselib.load_file(pathname, filename);
        o1 = data.Single_Shot_Occupations_SB1;
        o2 = data.Single_Shot_Occupations_SB2;
        o3 = data.Single_Shot_Occupations_SB3;
        o4 = data.Single_Shot_Occupations_SB4;
        
        X4 = o1 & o2 & o3 & o4;
        X3 = (o1 & o2 & o3) | (o2 & o3 & o4) | (o3 & o4 & o1) | (o4 & o1 & o2);
        for k = find(X4')'
            trial = fix(k/n_reps) + 1;
            rep = mod(k, n_reps);
            
            X3_trig = X3(trial,:);
            X3_trig = X3_trig( max(1, rep-50):min(n_reps,rep+50) );
            
            if sum(X3_trig) >= 3
                X4_event_date{end+1} = d;
                X4_event_file(end+1) = j;
                X4_event_trial(end+1) = trial;
                X4_event_rep(end+1) = rep;
                n_4X_events = n_4X_events + 1;

%                 autocorr_struct.q1 = add_to_autocorr(o1, autocorr_struct.q1, trial, rep, buffer_radius);
%                 autocorr_struct.q2 = add_to_autocorr(o2, autocorr_struct.q2, trial, rep, buffer_radius);
%                 autocorr_struct.q3 = add_to_autocorr(o3, autocorr_struct.q3, trial, rep, buffer_radius);
%                 autocorr_struct.q4 = add_to_autocorr(o4, autocorr_struct.q4, trial, rep, buffer_radius);
%                 autocorr_struct.q12 = add_to_autocorr(o1&o2, autocorr_struct.q12, trial, rep, buffer_radius);
%                 autocorr_struct.q34 = add_to_autocorr(o3&o4, autocorr_struct.q34, trial, rep, buffer_radius);
%                 autocorr_struct.q123 = add_to_autocorr(o1&o2&o3, autocorr_struct.q123, trial, rep, buffer_radius);
%                 autocorr_struct.q234 = add_to_autocorr(o2&o3&o4, autocorr_struct.q234, trial, rep, buffer_radius);
%                 autocorr_struct.q341 = add_to_autocorr(o3&o4&o1, autocorr_struct.q341, trial, rep, buffer_radius);
%                 autocorr_struct.q412 = add_to_autocorr(o4&o1&o2, autocorr_struct.q412, trial, rep, buffer_radius);
% 
%                 trial = randi(n_trials);
%                 rep = randi(n_reps);
%                 autocorr_struct.q1_baseline = add_to_autocorr(o1, autocorr_struct.q1_baseline, trial, rep, buffer_radius);
%                 autocorr_struct.q2_baseline = add_to_autocorr(o2, autocorr_struct.q2_baseline, trial, rep, buffer_radius);
%                 autocorr_struct.q3_baseline = add_to_autocorr(o3, autocorr_struct.q3_baseline, trial, rep, buffer_radius);
%                 autocorr_struct.q4_baseline = add_to_autocorr(o4, autocorr_struct.q4_baseline, trial, rep, buffer_radius);
%                 autocorr_struct.q12_baseline = add_to_autocorr(o1&o2, autocorr_struct.q12_baseline, trial, rep, buffer_radius);
%                 autocorr_struct.q34_baseline = add_to_autocorr(o3&o4, autocorr_struct.q34_baseline, trial, rep, buffer_radius);
%                 autocorr_struct.q123_baseline = add_to_autocorr(o1&o2&o3, autocorr_struct.q123_baseline, trial, rep, buffer_radius);
%                 autocorr_struct.q234_baseline = add_to_autocorr(o2&o3&o4, autocorr_struct.q234_baseline, trial, rep, buffer_radius);
%                 autocorr_struct.q341_baseline = add_to_autocorr(o3&o4&o1, autocorr_struct.q341_baseline, trial, rep, buffer_radius);
%                 autocorr_struct.q412_baseline = add_to_autocorr(o4&o1&o2, autocorr_struct.q412_baseline, trial, rep, buffer_radius);
            end
        end
        
        if mean(o1(:)) > 0.033
            hist_struct.q1 = add_to_hists(o1, hist_struct.q1, chunk_sizes);
        end
        if mean(o2(:)) > 0.06
            hist_struct.q2 = add_to_hists(o2, hist_struct.q2, chunk_sizes);
        end
        if mean(o3(:)) > 0.03
            hist_struct.q3 = add_to_hists(o3, hist_struct.q3, chunk_sizes);
        end
        if mean(o4(:)) > 0.046
            hist_struct.q4 = add_to_hists(o4, hist_struct.q4, chunk_sizes);
        end
        if (mean(o1(:)) > 0.033) && (mean(o2(:)) > 0.06)
            hist_struct.q12 = add_to_hists(o1&o2, hist_struct.q12, chunk_sizes);
        end
        if (mean(o3(:)) > 0.03) && (mean(o4(:)) > 0.046)
            hist_struct.q34 = add_to_hists(o3&o4, hist_struct.q34, chunk_sizes);
        end
        hist_struct.q123 = add_to_hists(o1&o2&o3, hist_struct.q123, chunk_sizes);
        hist_struct.q234 = add_to_hists(o2&o3&o4, hist_struct.q234, chunk_sizes);
        hist_struct.q341 = add_to_hists(o3&o4&o1, hist_struct.q341, chunk_sizes);
        hist_struct.q412 = add_to_hists(o4&o1&o2, hist_struct.q412, chunk_sizes);
        hist_struct.q1234 = add_to_hists(o1&o2&o3&o4, hist_struct.q1234, chunk_sizes);
        
        % average excited state
        avg_1_state(:,nf) = mean([o1(:) o2(:) o3(:) o4(:)]);
        
    end
end

close(f);

figure; plot(avg_1_state')
xlabel('file'); ylabel('Average Excess 1 state');

figure; hold on;
x_vals = -buffer_radius:buffer_radius;
plot(x_vals, autocorr_struct.q1/n_4X_events, 'DisplayName', 'Q1')
plot(x_vals, autocorr_struct.q2/n_4X_events, 'DisplayName', 'Q2')
plot(x_vals, autocorr_struct.q3/n_4X_events, 'DisplayName', 'Q3')
plot(x_vals, autocorr_struct.q4/n_4X_events, 'DisplayName', 'Q4')
plot(x_vals, autocorr_struct.q12/n_4X_events, 'DisplayName', 'Q12')
plot(x_vals, autocorr_struct.q34/n_4X_events, 'DisplayName', 'Q34')
plot(x_vals, (autocorr_struct.q123+autocorr_struct.q234 ...
            +autocorr_struct.q341+autocorr_struct.q412)/4/n_4X_events, 'DisplayName', 'Q_xyz')
plot(x_vals, autocorr_struct.q1234/n_4X_events, 'DisplayName', 'Q1234')
xlabel('Lag [100us]'); ylabel('Autocorrelation');
plot(x_vals, autocorr_struct.q1_baseline/n_4X_events, 'DisplayName', 'Q1_baseline')
plot(x_vals, autocorr_struct.q2_baseline/n_4X_events, 'DisplayName', 'Q2_baseline')
plot(x_vals, autocorr_struct.q3_baseline/n_4X_events, 'DisplayName', 'Q3_baseline')
plot(x_vals, autocorr_struct.q4_baseline/n_4X_events, 'DisplayName', 'Q4_baseline')
plot(x_vals, autocorr_struct.q12_baseline/n_4X_events, 'DisplayName', 'Q12_baseline')
plot(x_vals, autocorr_struct.q34_baseline/n_4X_events, 'DisplayName', 'Q34_baseline')
plot(x_vals, (autocorr_struct.q123_baseline+autocorr_struct.q234_baseline ...
            +autocorr_struct.q341_baseline+autocorr_struct.q412_baseline)/4/n_4X_events, 'DisplayName', 'Q_xyz_baseline')
plot(x_vals, autocorr_struct.q1234_baseline/n_4X_events, 'DisplayName', 'Q1234_baseline')
xlabel('Lag [100us]'); ylabel('Autocorrelation');

for n_reps_per_chunk = chunk_sizes
    rNum = ['r' num2str(n_reps_per_chunk)];
    figure; hold on;
    plot_poiss_hist(bin_edges, hist_struct.q1.(rNum), 'Q1')
    plot_poiss_hist(bin_edges, hist_struct.q2.(rNum), 'Q2')
    plot_poiss_hist(bin_edges, hist_struct.q3.(rNum), 'Q3')
    plot_poiss_hist(bin_edges, hist_struct.q4.(rNum), 'Q4')
    plot_poiss_hist(bin_edges, hist_struct.q12.(rNum), 'Q12')
    plot_poiss_hist(bin_edges, hist_struct.q34.(rNum), 'Q34')
    X3hist = hist_struct.q123.(rNum)+hist_struct.q234.(rNum)+hist_struct.q341.(rNum)+hist_struct.q412.(rNum);
    plot_poiss_hist(bin_edges, X3hist, 'Q_xyz')
    plot_poiss_hist(bin_edges, hist_struct.q1234.(rNum), 'Q1234')
    xlabel(['1 states in ' num2str(n_reps_per_chunk/10) 'ms chunk']); ylabel('Counts');
    set(gca, 'YScale', 'log');
    ylim([10^-1 inf]);
end

figure; hold on;
for n_reps_per_chunk = chunk_sizes
    rNum = ['r' num2str(n_reps_per_chunk)];
    plot_poiss_hist(bin_edges, hist_struct.q341.(rNum), ['Q341 ' num2str(n_reps_per_chunk)])
    xlabel(['1 states in chunk']); ylabel('Counts');
    set(gca, 'YScale', 'log');
    ylim([10^-1 inf]);
end


function [] = plot_poiss_hist(bin_edges, hist_bins, label)

histogram('BinEdges', bin_edges, 'BinCounts', hist_bins, 'DisplayName', label);
x = 0:length(hist_bins)-1;
errorbar(x, hist_bins, sqrt(hist_bins), '.');
n = sum(hist_bins);
avg = sum(x.*hist_bins)/n;
y = n*poisspdf(x,avg);
bins_with_counts = find(hist_bins);
max_bin = bins_with_counts(end);
deviation = ((hist_bins-y).^2)./y;
deviation = deviation(isfinite(deviation));
error_str = num2str( sum(deviation) / (max_bin + 1) );
plot(x, y, 'LineWidth', 2, 'DisplayName', [label ' fit ' error_str])

end

function [s] = add_to_hists(o, s, chunk_sizes)

bin_edges = -0.5:400.5;

for n_reps_per_chunk = chunk_sizes
    jumps_per_chunk = sum(reshape(o(:), n_reps_per_chunk, length(o(:))/n_reps_per_chunk));
    s.(['r' num2str(n_reps_per_chunk)]) = s.(['r' num2str(n_reps_per_chunk)]) ...
                                             + histcounts(jumps_per_chunk, bin_edges);
end

end


function [s] = add_to_autocorr(o, s, trial, rep, buffer_radius)

o_ac = o(trial,:);
o_ac = o_ac( max(1,rep-buffer_radius):min(length(o_ac),rep+buffer_radius) );
s = s + noiselib.crosscorrelate(o_ac,o_ac,buffer_radius);

end

