function data = maximizeIQContrast(data)
%maximizeIQContrast Find a rotated frame that maximize the data variations
%in the IQ-space.

for data_index = 1:length(data.dep)
    I_name = data.dep{data_index};
    
    if ~isempty(strfind(I_name, '_Std_Dev')) || isempty(strfind(I_name, 'I'))
        continue
    else
        Q_name = strrep(I_name, 'I', 'Q');
        if ~isfield(data, Q_name)
            continue
        end
        dep_name = strrep(I_name, 'I', 'Phase_Space_Shift');
        res_name = strrep(I_name, 'I', 'Residual');
    end

    I_rels = data.rels.(I_name);
    Q_rels = data.rels.(Q_name);
    
    if isempty(I_rels) || isempty(Q_rels) ||...
            length(I_rels) ~= length(Q_rels)
        continue
    else
        for k = 1:length(I_rels)
            if ~strcmp(I_rels{k}, Q_rels{k})
                error(['Dependecies do not match for ''', I_name, '''and''',...
                    Q_name, ''' data.'])
            end
        end
    end
    
    I = data.(I_name);
    Q = data.(Q_name);
    
    I = I - trimmean(I(:), 66);
    Q = Q - trimmean(Q(:), 66);
    phi = findAngle(I(:), Q(:));
    dep_vals = I * cos(phi) - Q * sin(phi);
    res_vals = I * sin(phi) + Q * cos(phi);
    
    data.(dep_name) = dep_vals;
    data.(res_name) = res_vals;
    data.rels.(dep_name) = I_rels;
    data.rels.(res_name) = I_rels;
    data.units.(dep_name) = data.units.(I_name);
    data.units.(res_name) = data.units.(I_name);
    data.dep{length(data.dep)+1} = dep_name;
    data.dep{length(data.dep)+1} = res_name; 
    
    if isfield(data, 'error') && isfield(data.error, I_name) &&...
        isfield(data.error, Q_name)
        data.error.(dep_name) = ...
                sqrt(data.error.(I_name).^2 * cos(phi)^2 +...
                     data.error.(Q_name).^2 * sin(phi)^2);
        data.error.(res_name) = ...
                sqrt(data.error.(I_name).^2 * sin(phi)^2 +...
                     data.error.(Q_name).^2 * cos(phi)^2);
    end
end
end

function phi = findAngle(I, Q)

phi = fminsearch(@costFun, 0);

    function SumSquaredError = costFun(phi)
        SumSquaredError = sum((I * sin(phi) + Q * cos(phi)).^2);
    end
end