function [tau_p, err_p, tau_r, err_r] = extractTimeConstants(t, n_qp)

% Extract poisoning time constant.
t_p = t(t < 0) - min(t);
n_p = n_qp(t < 0);

f_p = fit(t_p(:), n_p(:),...
        'a * (1 - exp(-b * x))', 'StartPoint', [max(n_p), 1]);
    
tau_p = 1 / f_p.b;
ci_p = confint(f_p);
err_p = max([abs(1 / ci_p(1, 2) - 1 / f_p.b),...
             abs(1 / ci_p(2, 2) - 1 / f_p.b)]);

% Extract recombination time constant.
t_r = t(t >= 0);
n_r = t(t >= 0);

f_r = fit(t_r(:), n_r(:),...
        'a * exp(-b * x)', 'StartPoint', [max(n_r), 1]);
    
tau_r = 1 / f_r.b;
ci_r = confint(f_r);
err_r = max([abs(1 / ci_r(1, 2) - 1 / f_r.b),...
             abs(1 / ci_r(2, 2) - 1 / f_r.b)]);
end