function data = importIVCharacteristic(filename)

data = [];

if nargin ~= 1
    disp(['Error using ', mfilename('fullpath'), ': ', 'Function requires one input parameter'])
    return
end

[fid, msg] = fopen(filename, 'r');

if fid == -1
    disp(['Error using ', mfilename('fullpath'), ': ', msg]);
    return
end

line = fgetl(fid);
if ~ischar(line)
    disp(['Error using ', mfilename('fullpath'), ': ', 'Datafile ', filename, ' is empty']);
    return
end

% data.Current_Bias_Resistance = 1000; % Ohm (unless specified in the header)
% data.Voltage_Gain = 1000; % dimensionless (unless secified in the header)

for k = 1:sscanf(line, '%d', 1)-1
    line = fgetl(fid);
    if ~ischar(line)
        disp(['Error using ', mfilename('fullpath'), ': ', 'Datafile ', filename, ' is too short']);
        return
    end
    data.comment{k} = line;
%     pos = strfind(line, ':');
%     if ~isempty(pos) && length(pos) == 1
%         if ~isempty(strfind(line, 'Resistance')) || ~isempty(strfind(line, 'Ohm'))
%             data.Current_Bias_Resistance = sscanf(line(pos+1:end), '%f', 1);
%         end
%         if ~isempty(strfind(line, 'Gain'))
%             data.Voltage_Gain = sscanf(line(pos+1:end), '%f', 1);
%         end
%     end
end

val = (fscanf(fid, '%f', [2, inf]))';
fclose(fid);

data.Voltage = val(:, 1);
data.Current = val(:, 2);