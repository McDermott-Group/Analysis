function [unwrap_es] = unwrap_charge(vs, wrap_voltage, dedv, unwrap_mult)

unwrap_vs = vs;
wrap_value = 0;
for i = 1:length(vs)
    unwrap_vs(i) = vs(i) + unwrap_mult*wrap_value/dedv;
    if vs(i) > wrap_voltage% + 7.5
        wrap_value = wrap_value + 1;
    end
    if vs(i) < -wrap_voltage% + 7.5
        wrap_value = wrap_value - 1;
    end
end
unwrap_es = unwrap_vs * dedv;

end