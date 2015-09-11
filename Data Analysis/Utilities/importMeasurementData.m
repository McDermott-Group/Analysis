function data = importMeasurementData(filename)
%importMeasurementData Import data from a text file created with the JPM
%Experiment code.
%
%   DATA = importMeasurementData(FILENAME) reads data from a file specified
%   by FILENAME and returns structure DATA containing the data from
%   the FILENAME file. Only 1D and 2D dataseta are supported.

if nargin ~= 1
    error('Function requires only one input parameter: FILENAME.');
end

[fid, msg] = fopen(filename, 'r');

if fid == -1
    error(msg);
end

data.Filename = filename;

line = fgetl(fid);
if ~ischar(line)
    error(['Datafile ', filename, ' is empty. The first line should be the name of the experiment.']);
else
    data.Experiment_Name = line;
end

line = fgetl(fid);
if ~ischar(line)
    error(['Datafile ', filename, ' is empty. The second line should be the timestamp.']);
else
    data.Timestamp = line;
end

line = fgetl(fid);
if isempty(strfind(line, '====Experiment Parameters===='))
    error(['"====Experiment Parameters====" line is expected in ', filename, '.']);
end

comment_counter = 0;
line = fgetl(fid);
while isempty(strfind(line, '====Sweep Variables===='))
    if ~ischar(line)
        error(['"====Sweep Variables====" line is expected in ' , filename, '.']);
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
                data.Comments{comment_counter} = line(4:end);
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
        error(['"====Data Variables====" line is expected in ' , filename, '.']);
    end
    if ~isempty(line)
        pos = strfind(line, ':');
        while length(pos) ~= 1
            if ~ischar(line)
                error(['Cannot recognize a sweep variable name in ', filename, '.']);
            end
            line = fgetl(fid);
            pos = strfind(line, ':');
        end
        ubra = strfind(line, '[');
        uket = strfind(line, ']');
        if isempty(ubra) || isempty(uket)
            error(['Sweep variables are not properly specified in ', filename, '.']);
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
                error(['Data variable dependencies are not properly specified in ', filename, '.']);
            else
                relationships = cell(1, length(rels_pos)/2);
                for k = 1:2:length(rels_pos)
                    relationships{(k+1)/2} = strrep(deps_line(rels_pos(k)+1:rels_pos(k+1)-1), ' ', '_');
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
            error('Data variables that have more than four dimensions are not yet supported.')
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

fclose(fid);