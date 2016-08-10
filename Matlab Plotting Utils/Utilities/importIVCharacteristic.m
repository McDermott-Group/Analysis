function data = importIVCharacteristic(filename)

data = [];

if nargin ~= 1
    disp(['Error using ', mfilename('fullpath'), ': ',...
        'Function requires one input parameter'])
    return
end

[fid, msg] = fopen(filename, 'r');

if fid == -1
    disp(['Error using ', mfilename('fullpath'), ': ', msg]);
    return
end

line = fgetl(fid);
if ~ischar(line)
    disp(['Error using ', mfilename('fullpath'), ': ', 'Datafile ',...
        filename, ' is empty']);
    return
end

for k = 1:sscanf(line, '%d', 1)-1
    line = fgetl(fid);
    if ~ischar(line)
        disp(['Error using ', mfilename('fullpath'), ': ', 'Datafile ',...
            filename, ' is too short.']);
        return
    end
    data.comment{k} = line;
end

val = (fscanf(fid, '%f', [2, inf]))';
fclose(fid);

data.Voltage = val(:, 1);
data.Current = val(:, 2);