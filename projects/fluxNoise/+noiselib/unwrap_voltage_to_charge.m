function [unwrap_es] = unwrap_voltage_to_charge(vs, wrap_voltage, dedv)

unwrap_vs = vs;
wrap_value = 0;
for i = 1:length(vs)
    unwrap_vs(i) = vs(i) + wrap_value/dedv;
    if vs(i) > wrap_voltage% + 7.5
        wrap_value = wrap_value + 1;
    end
    if vs(i) < -wrap_voltage% + 7.5
        wrap_value = wrap_value - 1;
    end
end
unwrap_es = unwrap_vs * dedv;

% take into account aliasing
delta_q = (unwrap_es(2:end) - unwrap_es(1:end-1));
delta_q = noiselib.alias(delta_q, 0.5);
first_point = noiselib.alias(unwrap_es(1), 0.5);
unwrap_es = [[first_point]; first_point+cumsum(delta_q)];

end

