function [tau_p, err_p, tau_r, err_r] = extractTimeConstants(t, n_qp, plot_flag)

if ~exist('plot_flag', 'var')
    plot_flag = false;
end

% Extract poisoning time constant.
t_p = t(t < 0) - min(t);
n_p = n_qp(t < 0);

[~, idx] = max(n_qp);
t_start = t_p(idx);

f_p = fit(t_p(:), n_p(:),...
        'a * (1 - exp(-b * x))', 'StartPoint', [max(n_p), 10 / t_start]);

tau_p = 1 / f_p.b;
ci_p = confint(f_p);
err_p = max([abs(1 / ci_p(1, 2) - 1 / f_p.b),...
             abs(1 / ci_p(2, 2) - 1 / f_p.b)]);

if plot_flag
    figure
    hold on
    plot(t_p, n_p, 'LineWidth', 2)
    plot(t_p, f_p.a * (1 - exp(-t_p / tau_p)), 'LineWidth', 2)
    xlabel('Time (\tau_0)', 'FontSize', 14)
    ylabel('n_{\rm qp} / n_{\rm cp}', 'FontSize', 14)
    title('Poisoning Fit')
    legend('simulation', 'fit')
    grid on
    grid minor
    axis tight
end

% Extract recombination time constant.
t_r = t(t > 0);
n_r = n_qp(t > 0);

f_r = fit(t_r(:), n_r(:),...
        'a * exp(-b * x)', 'StartPoint', [max(n_r), f_p.b]);
    
tau_r = 1 / f_r.b;
ci_r = confint(f_r);
err_r = max([abs(1 / ci_r(1, 2) - 1 / f_r.b),...
             abs(1 / ci_r(2, 2) - 1 / f_r.b)]);

if plot_flag
    figure
    hold on
    plot(t_r, n_r, 'LineWidth', 2)
    plot(t_r, f_r.a * exp(-t_r / tau_r), 'LineWidth', 2)
    xlabel('Time (\tau_0)', 'FontSize', 14)
    ylabel('n_{\rm qp} / n_{\rm cp}', 'FontSize', 14)
    title('Relaxation Fit')
    legend('simulation', 'fit')
    grid on
    grid minor
    axis tight
end

end