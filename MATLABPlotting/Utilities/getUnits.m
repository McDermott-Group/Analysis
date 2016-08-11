function units = getUnits(data, var_name)
%getUnits   Return variable units.
%
%   UNITS = getUnits(DATA, VAR_NAME) returns a string that specify the
%   variable units. The units are placed in brackets. If no units are
%   specified or the variable is unitless, an empty string is returned.

if isfield(data.units, var_name) && ~isempty(data.units.(var_name))
    units = [' (', data.units.(var_name), ')'];
else
    units = '';
end

end

