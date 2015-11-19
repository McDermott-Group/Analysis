function createFigure(position)
%createFigure   Create a figure.
%
%   createFigure(POSITION) creates a figure in the POSITION that is defined
%   relatively to the screen size.

if exist('position', 'var')
    if ischar(position) && strcmpi(position, 'right')
        str_position = 'right';
    elseif isnumeric(position) && length(position(:)) == 4
        str_position = 'custom';
    end
end
if ~exist('str_position', 'var')
    str_position = 'left';
end

if strcmp(str_position, 'left') 
    position = [.10, .40, .75, .50];
elseif strcmp(str_position, 'right') 
    position = [.85, .40, .75, .50];
end

scrsz = get(0, 'ScreenSize');
figure('Position', position(:) * scrsz(4));

end

