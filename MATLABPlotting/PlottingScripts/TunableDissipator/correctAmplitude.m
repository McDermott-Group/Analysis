function correctAmplitude
%correctAmplitude Correct the amplitude by shifting, rescaling and
%offsetting the in-phase and quadrature components.

% Select a file.
data = loadMeasurementData;
if isempty(fields(data))
    return
end

for data_index = 1:length(data.dep)
    I_name = data.dep{data_index};
    
    if ~isempty(strfind(I_name, '_Std_Dev')) ||...
            isempty(strfind(I_name, 'I')) ||...
            length(data.(I_name)) < 2 ||...
            length(data.rels.(I_name)) ~= 1
        continue
    else
        Q_name = strrep(I_name, 'I', 'Q');
        if ~isfield(data, Q_name)
            continue
        end
        dep_name = 'Corrected_Amplitude';
    end

    I = data.(I_name);
    Q = data.(Q_name);

    fitted_params = findParams(I(:), Q(:));
    
    dep_vals = correctedAmplitude(I, Q, fitted_params);

    data.(dep_name) = dep_vals;
    data.rels.(dep_name) = data.rels.(I_name);
    data.units.(dep_name) = data.units.(I_name);
    data.dep{length(data.dep)+1} = dep_name;

    plotDataVar(data, dep_name);
    hold on
    plot(data.(data.rels.(I_name){1}), hypot(I, Q))
    hold off
    legend({'corrected', 'original'}) 
    
    [~, filename, ~] = fileparts(data.Filename);
    saveMeasData(data, [filename, '_corrected_amplitude'])
end
end

function fitted_params = findParams(I, Q)

fitted_params = fminsearch(@costFun, [mean(I), mean(Q), 1, 0]);

    function SumSquaredError = costFun(params)
        A = correctedAmplitude(I, Q, params);
        SumSquaredError = sum(hypot(1, diff(A)));
    end
end

function A = correctedAmplitude(I, Q, x)
    I0 = x(1);
    Q0 = x(2);
    ratio = x(3);
    shift = x(4);

    t = 1:length(I);
    FQ = griddedInterpolant(t, Q, 'spline');
    Q_interpolated = FQ(t + shift);
    A = hypot(I - I0,...
             ratio * (Q_interpolated(:) - Q0));
end