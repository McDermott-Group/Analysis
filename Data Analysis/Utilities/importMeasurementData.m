function data = importMeasurementData(filename)
%importMeasurementData Import data from a text file created with the JPM
%Experiment code.
%
%   DATA = importMeasurementData(FILENAME) reads data from a file specified
%   by FILENAME and returns structure DATA containing the data from
%   the FILENAME file. Only 1D and 2D dataseta are supported.

    if nargin ~= 1
        error('Function requires one and only one input parameter: FILENAME.');
    end

    [~, ~, ext] = fileparts(filename);
    if strcmp(ext, '.mat')
        data = importMat_v0p0(filename);
    elseif strcmp(ext, '.txt')
        [fid, msg] = fopen(filename, 'r');

        if fid == -1
            error(msg);
        end

        line = fgetl(fid);
        if ~ischar(line)
            error(['Datafile ', filename, ' is empty. The first line should be',...
                'either the format version (or the name of the experiment',...
                ' for the  data sets).']);
        elseif strcmp(line, 'Format Version: 0.1')
            data = importTxt_v0p1(filename, fid, line);
        else
            data = importTxt_v0p0(filename, fid, line);
        end
        fclose(fid);
    else
        error(['Cannot read ', ext(2:end), '-files.'])
    end
end

function data = importMat_v0p0(filename)
    [~, fn, ~] = fileparts(filename);
    data = load(filename);
    data = data.(fn);
    
    data.Filename = filename;
    data.Timestamp = data.Time;
    data = rmfield(data, 'Time');
    
    data.units = data.Units;
    data = rmfield(data, 'Units');
    
    if isfield(data, 'Distr')
        data.distr = data.Distr;
        data = rmfield(data, 'Distr');
    else
        data.distr = struct();
    end
    
    data.rels = data.Depend;
    data = rmfield(data, 'Depend');
    
    comments = data.Comments;
    data.Comments = cell(1);
    data.Comments{1} = comments; 
    
    rels = fieldnames(data.rels);
    for q = 1:length(rels)
        rels_str = data.rels.(rels{q});
        rels_pos = strfind(rels_str, char(39));    % char(39) is a single quotation mark
        if mod(length(rels_pos), 2) ~= 0
            error(['Data variable dependencies are not properly',...
                'specified in ', filename, '.']);
        end
        relationships = cell(1, length(rels_pos)/2);
        for k = 1:2:length(rels_pos)
            relationships{(k+1)/2} = strrep(rels_str(rels_pos(k)+1:rels_pos(k+1)-1),...
                ' ', '_');
        end
        data.rels.(rels{q}) = relationships;
    end
            
    expt_vars = fieldnames(data.ExptVars);
    for k = 1:length(expt_vars)
        data.(expt_vars{k}) = data.ExptVars.(expt_vars{k});
    end
    data = rmfield(data, 'ExptVars');
    
    dep_counter = 0;
    indep_counter = 0;
    rels = fieldnames(data.rels);
    vars = fieldnames(data.Data);
    for k = 1:length(vars)
        is_dep = true;
        for q = 1:length(rels)
            for p = 1:length(data.rels.(rels{q}))
                if strcmp(data.rels.(rels{q}){p}, vars{k})
                    is_dep = false;
                    break
                end
            end
        end
        if is_dep
            dep_counter = dep_counter + 1;
            data.dep{dep_counter} = vars{k};
        else
            indep_counter = indep_counter + 1;
            data.indep{indep_counter} = vars{k};
        end
        data.(vars{k}) = data.Data.(vars{k});
    end
    data = rmfield(data, 'Data');

    for k = 1:length(data.dep)
        if ~isfield(data.distr, data.dep{k})
            data.distr.(data.dep{k}) = '';
        end
    end
end

function data = importTxt_v0p1(filename, fid, first_line)
    data.Filename = filename;
    data.Text_Format_Version = first_line(strfind(first_line, ': ')+2:end);
    
    line = fgetl(fid);
    if ~ischar(line)
        error(['Datafile ', filename, ' is empty. The first line should',...
            ' be the experiment name.']);
    else
        data.Experiment_Name = line;
    end
    
    line = fgetl(fid);
    if ~ischar(line)
        error(['Datafile ', filename, ' is empty. The second line should',...
            ' be the timestamp.']);
    else
        data.Timestamp = line;
    end

    line = fgetl(fid);
    if isempty(strfind(line, '====Experiment Parameters===='))
        error(['"====Experiment Parameters====" line is expected in ',...
            filename, '.']);
    end

    comment_counter = 0;
    line = fgetl(fid);
    while isempty(strfind(line, '====Sweep Variables===='))
        if ~ischar(line)
            error(['"====Sweep Variables====" line is expected in ',...
                filename, '.']);
        end
        if ~isempty(line)
            pos = strfind(line, ':');
            if isempty(pos) || length(pos) > 1
                comment_counter = comment_counter + 1;
                if length(pos) > 1
                    data.Comments{comment_counter} = line(pos(1)+1:end);
                else
                    data.Comments{comment_counter} = line;
                end
            else
                fieldname = strrep(line(1:pos-1), ' ', '_');
                line = line(pos+1:end);
                if strcmp(fieldname, 'Comments')
                    comment_counter = comment_counter + 1;
                    data.Comments{comment_counter} = line(2:end);
                else
                    data.(fieldname) = sscanf(line, '%f %*s');
                end
                ubra = strfind(line, '[');
                uket = strfind(line, ']');
                if ~isempty(ubra) && ~isempty(uket)
                    data.units.(fieldname) = line(ubra+1:uket-1);
                elseif ~strcmp(fieldname, 'Comments')
                    data.units.(fieldname) = '';
                end
            end
        end
        line = fgetl(fid);
    end

    line = fgetl(fid);
    indep_counter = 0;
    while isempty(strfind(line, '====Data Variables===='))
        if ~ischar(line)
            error(['"====Data Variables====" line is expected in ',...
                filename, '.']);
        end
        if ~isempty(line)
            if ~isempty(strfind(line, 'Name: '))
                pos = strfind(line, char(39));
                if length(pos) ~= 2
                    error(['Cannot recognize an independent variable ',...
                        'name in ', filename, '.']);
                end
                indep_counter = indep_counter + 1;
                data.indep{indep_counter} = strrep(line(pos(1)+1:pos(2)-1),...
                    ' ', '_');
                data.units.(data.indep{indep_counter}) = '';
            elseif ~isempty(strfind(line, 'Units: '))
                pos = strfind(line, ': ');
                data.units.(data.indep{indep_counter}) = line(pos+2:end);
            elseif ~isempty(strfind(line, 'Type: ')) &&... 
                    isempty(strfind(line, 'independent'))
                error(['Unknown type detected among the sweep ',...
                    '(independent) variables.']);
            elseif ~isempty(strfind(line, 'Size: '))
                pos = strfind(line, ':');
                sz = sscanf(line(pos+1:end), '%*c%f');
                if isempty(sz)
                    error('Sweep variables should not be empty.')
                end
                if length(sz) > 1
                    error(['Sweep variables should not have more than ',...
                        'one dimension.'])
                end
                data.(data.indep{indep_counter}) = fscanf(fid, '%f', sz);
            end
        end
        line = fgetl(fid);
    end

    dep_counter = 0;
    while ischar(line)
        if ~isempty(strfind(line, 'Name:'))
            pos = strfind(line, char(39));
            if length(pos) ~= 2
                error(['Cannot recognize an dependent variable ',...
                    'name in ', filename, '.']);
            end
            dep_counter = dep_counter + 1;
            data.dep{dep_counter} = strrep(line(pos(1)+1:pos(2)-1),...
                ' ', '_');
            data.units.(data.dep{dep_counter}) = '';
            data.distr.(data.dep{dep_counter}) =  '';
        elseif ~isempty(strfind(line, 'Units: '))
            pos = strfind(line, ': ');
            data.units.(data.dep{dep_counter}) = line(pos+2:end);
        elseif ~isempty(strfind(line, 'Type: ')) &&...
                isempty(strfind(line, ' dependent'))
            error(['Unknown type detected among the data ',...
                '(dependent) variables.']);
        elseif ~isempty(strfind(line, 'Distribution:'))     
            pos = strfind(line, ': ');
            data.distr.(data.dep{dep_counter}) = line(pos+2:end);
        elseif ~isempty(strfind(line, 'Dependencies: '))
            pos = strfind(line, ': ');
            deps_line = line(pos+2:end);
            rels_pos = strfind(deps_line, char(39));    % char(39) is a single quotation mark
            if mod(length(rels_pos), 2) ~= 0
                error(['Data variable dependencies are not properly',...
                    'specified in ', filename, '.']);
            end
            relationships = cell(1, length(rels_pos)/2);
            for k = 1:2:length(rels_pos)
                relationships{(k+1)/2} = strrep(deps_line(rels_pos(k)+1:rels_pos(k+1)-1),...
                    ' ', '_');
            end
            data.rels.(data.dep{dep_counter}) = relationships;
        elseif ~isempty(strfind(line, 'Size: '))
            pos = strfind(line, ': ');
            sz = sscanf(line(pos+1:end), '%*c%f');
            if isempty(sz)
                error('Data variables should not be empty.')
            end
            if length(sz) > 4
                error(['Data variables that have more than four ',...
                    'dimensions are not yet supported.'])
            end
            if length(sz) == 1
                data.(data.dep{dep_counter}) = fscanf(fid, '%f', sz);
            elseif length(sz) == 2
                data.(data.dep{dep_counter}) = (fscanf(fid, '%f', [sz(2), sz(1)]))';
            elseif length(sz) == 3
                temp_data = nan(sz');
                for k = 1:sz(1)
                    temp_data(k, :, :) = (fscanf(fid, '%f', [sz(3), sz(2)]))';
                end
                data.(data.dep{dep_counter}) = temp_data;
            elseif length(sz) == 4
                temp_data = nan(sz');
                for k1 = 1:sz(1)
                    for k2 = 1:sz(2)
                        temp_data(k1, k2, :, :) = (fscanf(fid, '%f', [sz(4), sz(3)]))';
                    end
                end
                data.(data.dep{dep_counter}) = temp_data;
            end
        end
        line = fgetl(fid);
    end
end

function data = importTxt_v0p0(filename, fid, first_line)
    data.Filename = filename;
    data.Text_Format_Version = '0.0';
    data.Experiment_Name = first_line;
    
    line = fgetl(fid);
    if ~ischar(line)
        error(['Datafile ', filename, ' is empty. The second line ',...
            'should be the timestamp.']);
    else
        data.Timestamp = line;
    end

    line = fgetl(fid);
    if isempty(strfind(line, '====Experiment Parameters===='))
        error(['"====Experiment Parameters====" line is expected in ',...
            filename, '.']);
    end

    comment_counter = 0;
    line = fgetl(fid);
    while isempty(strfind(line, '====Sweep Variables===='))
        if ~ischar(line)
            error(['"====Sweep Variables====" line is expected in ',...
                filename, '.']);
        end
        if ~isempty(line)
            pos = strfind(line, ':');
            if isempty(pos) || length(pos) > 1
                comment_counter = comment_counter + 1;
                if length(pos) > 1
                    data.Comments{comment_counter} = line(pos(1)+1:end);
                else
                    data.Comments{comment_counter} = line;
                end
            else
                fieldname = strrep(line(1:pos-1), ' ', '_');
                line = line(pos+1:end);
                if strcmp(fieldname, 'Comments')
                    comment_counter = comment_counter + 1;
                    data.Comments{comment_counter} = line(2:end);
                else
                    data.(fieldname) = sscanf(line, '%f %*s');
                end
                ubra = strfind(line, '[');
                uket = strfind(line, ']');
                if ~isempty(ubra) && ~isempty(uket)
                    data.units.(fieldname) = line(ubra+1:uket-1);
                elseif ~strcmp(fieldname, 'Comments')
                    data.units.(fieldname) = '';
                end
            end
        end
        line = fgetl(fid);
    end

    line = fgetl(fid);
    indep_counter = 0;
    while isempty(strfind(line, '====Data Variables===='))
        if ~ischar(line)
            error(['"====Data Variables====" line is expected in ',...
                filename, '.']);
        end
        if ~isempty(line)
            pos = strfind(line, ':');
            while length(pos) ~= 1
                if ~ischar(line)
                    error(['Cannot recognize a sweep variable name in ',...
                        filename, '.']);
                end
                line = fgetl(fid);
                pos = strfind(line, ':');
            end
            ubra = strfind(line, '[');
            uket = strfind(line, ']');
            if isempty(ubra) || isempty(uket)
                error(['Sweep variables are not properly specified in ',...
                    filename, '.']);
            end
            indep_counter = indep_counter + 1;
            if ubra(1) < pos && uket(1) < pos
                data.indep{indep_counter} = strrep(line(1:ubra-2), ' ', '_');
                data.units.(data.indep{indep_counter}) = line(ubra+1:uket-1);
            else
                data.indep{indep_counter} = strrep(line(1:pos-1), ' ', '_');
                data.units.(data.indep{indep_counter}) = '';
            end
            sz = sscanf(line(pos+1:end), '%*c%f');
            if isempty(sz)
                error('Sweep variables should not be empty.')
            end
            if length(sz) > 1
                error('Sweep variables should not have more than one dimension.')
            end
            data.(data.indep{indep_counter}) = fscanf(fid, '%f', sz);
        end
        line = fgetl(fid);
    end

    dep_counter = 0;
    while ischar(line)      % At the first interation line contains '====Data Variables===='.
        line = fgetl(fid);
        if ~isempty(line)   % Check whether the end of file is reached.
            pos = strfind(line, ':');
            if isempty(pos) % Get another line.
                continue
            end
            ubra = strfind(line, '[');
            uket = strfind(line, ']');
            if isempty(ubra) || isempty(uket)
                error(['Data variables are not properly specified in ', filename, '.']);
            end
            dep_counter = dep_counter + 1;
            if ubra(1) < pos(1) && uket(1) < pos(1)
                data.dep{dep_counter} = strrep(line(1:ubra-2), ' ', '_');
                data.units.(data.dep{dep_counter}) = line(ubra+1:uket-1);
            else
                data.dep{dep_counter} = strrep(line(1:pos(1)-1), ' ', '_');
                data.units.(data.dep{dep_counter}) = '';
            end

            rels_pos = strfind(line, '::');
            relationships = {};
            if ~isempty(rels_pos)
                sz = sscanf(line(pos(1)+1:rels_pos-1), '%*c%f');
                deps_line = line(rels_pos+2:end);
                rels_pos = strfind(deps_line, char(39));     % single quote
                if mod(length(rels_pos), 2) ~= 0
                    error(['Data variable dependencies are not properly specified in ',...
                        filename, '.']);
                else
                    relationships = cell(1, length(rels_pos)/2);
                    for k = 1:2:length(rels_pos)
                        relationships{(k+1)/2} = strrep(deps_line(rels_pos(k)+1:rels_pos(k+1)-1),...
                        ' ', '_');
                    end
                end
            else
                sz = sscanf(line(pos+1:end), '%*c%f');
            end
            data.rels.(data.dep{dep_counter}) = relationships;

            distr_pos = strfind(line, ':::');
            if ~isempty(distr_pos)
                data.distr.(data.dep{dep_counter}) = line(distr_pos+4:end-1);
            else
                data.distr.(data.dep{dep_counter}) =  '';
            end

            if isempty(sz)
                error('Data variables should not be empty.')
            end
            if length(sz) > 4
                error(['Data variables that have more than four ',...
                    'dimensions are not yet supported.'])
            end
            if length(sz) == 1
                data.(data.dep{dep_counter}) = fscanf(fid, '%f', sz);
            elseif length(sz) == 2
                data.(data.dep{dep_counter}) = (fscanf(fid, '%f', [sz(2), sz(1)]))';
            elseif length(sz) == 3
                temp_data = nan(sz');
                for k = 1:sz(1)
                    temp_data(k, :, :) = (fscanf(fid, '%f', [sz(3), sz(2)]))';
                end
                data.(data.dep{dep_counter}) = temp_data;
            elseif length(sz) == 4
                temp_data = nan(sz');
                for k1 = 1:sz(1)
                    for k2 = 1:sz(2)
                        temp_data(k1, k2, :, :) = (fscanf(fid, '%f', [sz(4), sz(3)]))';
                    end
                end
                data.(data.dep{dep_counter}) = temp_data;
            end
        end
    end
end