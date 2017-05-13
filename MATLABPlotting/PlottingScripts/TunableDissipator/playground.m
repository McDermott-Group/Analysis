function playground
%playground Test of some stupid ideas.

data = loadMeasurementData;

I_name = 'It';
Q_name = 'Qt';
dep_name = 'Corrected_Amplitude';

I = data.(I_name);
Q = data.(Q_name);

params = findParams(I(:), Q(:))

t = 1:length(I);
FQ = griddedInterpolant(t, Q);
Q_interpolated = FQ(t + params(4));

dep_vals = hypot(I - params(1), params(3) * (Q_interpolated(:) - params(2)));

data.(dep_name) = dep_vals;
data.rels.(dep_name) = data.rels.(I_name);
data.units.(dep_name) = data.units.(I_name);
data.dep{length(data.dep)+1} = dep_name;

plotDataVar(data, dep_name);

end

function x = findParams(I, Q)

x = fminsearch(@costFun, [mean(I), mean(Q), 1, 0])

    function SumSquaredError = costFun(x)

        I0 = x(1);
        Q0 = x(2);
        ratio = x(3);
        shift = x(4);
        
        t = 1:length(I);
        FQ = griddedInterpolant(t, Q);
        Q_interpolated = FQ(t + shift);
        A = hypot(I - I0, ratio * (Q_interpolated(:) - Q0));
        SumSquaredError = sum(hypot(1, diff(A)));
    end
end