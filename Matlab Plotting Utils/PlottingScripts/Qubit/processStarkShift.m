function processStarkShift(data, chi_0, alpha, f_cavity)
%processStarkShift(DATA) Fit photon number to DAC amplitude assuming known
% anharmonicity of qubit and assuming a three-level system.
%   processStarkShift(DATA, CHI_0, ALPHA, F_CAVITY) fits photon number to
%   DAC Stark amplitude specified in the data structure DATA, assuming
%   CHI_0 shift, given in MHz, and anharmonicity of the qubit ALPHA,
%   specified in MHz (the value should be negative). F_CAVITY is the bare
%   cavity bare frequency.
%
%   For details, see M. Boissonneault et al., PRL 105, 100504 (2010):
%   http://journals.aps.org/prl/pdf/10.1103/PhysRevLett.105.100504

stark_amplitude_name = 'Stark_Amplitude';
if isfield(data, stark_amplitude_name)
    stark_amplitude = data.(stark_amplitude_name);
else
    error('Stark Amplitude data is not found.')
end

stark_frequency_name = 'Extracted_Resonance_Frequency';
if isfield(data, stark_frequency_name)
    fitted_stark_frequency = data.(stark_frequency_name);
else
    stark_frequency_name = 'Stark_Frequency';
    if isfield(data, stark_frequency_name)
        fitted_stark_frequency = data.(stark_frequency_name);
    else
        error(['Neither Extracted Resonance Freqeuency nor ',...
            'Stark Frequency data is not found.'])
    end
end
unit = data.units.(stark_frequency_name);
if strcmp(unit, 'GHz')
    fitted_stark_frequency = fitted_stark_frequency * 1e9;
else
    error(['Cannot convert ', unit, ' to Hz.'])
end

if exist('alpha', 'var')
    alpha = alpha * 1e6; % Hz, quibit anharmonicity
else
    error('Qubit anharmonicity is not specified.')
end

if exist('chi_0', 'var')
    chi_0 = chi_0 * 1e6; % Hz, low power chi shift
else
    error('Chi shift is not specified.')
end

if exist('f_cavity', 'var')
    f_cavity = f_cavity * 1e9; % Hz, bare cavity frequency
else
    error('Bare cavity frequency is not specified.')
end

w_10 = 2 * pi * fitted_stark_frequency(1); % s^-1, qubit 10 frequency
w_r  = 2 * pi * f_cavity; % s^-1
w_21 = 2 * pi * (fitted_stark_frequency(1) + alpha); % s^-1

% Detunings.
delta_0 = w_10 - w_r;
delta_1 = w_21 - w_r;

g_10 = sqrt(2 * pi * chi_0 * abs(delta_0));
g_21 = sqrt(2) * g_10;

lambda_0 = -g_10 / delta_0;
lambda_1 = -g_21 / delta_1;

Chi_0 = g_10^2 / delta_0;
Chi_1 = g_21^2 / delta_1;

g2_0 = lambda_0 * lambda_1 * (delta_1 - delta_0);
lambda2_0 = -g2_0 / (delta_1 + delta_0);

% Stark and kerr coefficients.
S_0 = -Chi_0 + (-Chi_0 * lambda_1^2 + 3 * Chi_1 * lambda_0^2) / 4 -...
    g2_0 * lambda2_0;
S_1 = Chi_0 * (1 - lambda_1^2) - Chi_1 * (1 - lambda_0^2) -...
    2 * Chi_0 * lambda_0^2;

K_0 = (Chi_0 * lambda_1^2 - 3 * Chi_1 * lambda_0^2) / 4 -...
    g2_0 * lambda2_0;
K_1 = (Chi_1 - Chi_0) * (lambda_1^2 + lambda_0^2);

n_crit = 1 /(4 * lambda_0^2);

shift = 2 * pi * fitted_stark_frequency - w_10;

    function sse = stark_fit(x)
        fit_curve = (S_1 - S_0) * x   * stark_amplitude.^2 +...
                    (K_1 - K_0) * x^2 * stark_amplitude.^4;
        err = fit_curve - shift;
        sse = sum(err.^2);
    end

stark_factor = fminsearch(@stark_fit, 0.1);

SSA = linspace(min(stark_amplitude), max(stark_amplitude), 1001);
SS_fit = (S_1 - S_0) * stark_factor * SSA.^2 +...
         (K_1 - K_0) * stark_factor * SSA.^4;
n_fit = stark_factor * SSA.^2;

createFigure([.2, .1, .75, .75]);
subplot(2, 1, 1)
[~, filename, ext] = fileparts(data.Filename);
title({[strrep(filename, '_', '\_'), ext, ' [', data.Timestamp, ']'],...
    ['Fit to Stark Shift Data: number of photons = ',...
    num2str(stark_factor, 4), ' x (DAC Amplitude)^2'],...
    ['n_{crit} = ', num2str(n_crit, 4), ' photons; ',...
    'Qubit Anharmonicity = ', num2str(alpha / 1e6), ' MHz; ',...
    '\chi_0 = ', num2str(chi_0/1e6), ' MHz']}, 'FontSize', 10)
hold on
plot(stark_amplitude, shift / pi / 2e6, 'r.', 'MarkerSize', 15)
plot(SSA, SS_fit / pi / 2e6, 'b-', 'LineWidth', 2)
hold off
ylabel('Stark Shift (MHz)', 'FontSize', 14)
grid on

subplot(2, 1, 2)
plot(SSA, n_fit, 'b-', 'LineWidth', 2)
ylabel('Number of Photons', 'FontSize', 14)
xlabel([strrep(stark_amplitude_name, '_', ' '),...
    getUnits(data, stark_amplitude_name)], 'FontSize', 14)
grid on

end