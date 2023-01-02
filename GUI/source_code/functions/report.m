function report(message, handles, col)
if nargin < 3; col = 'b'; end
    set(handles.edit7, 'String', message, 'ForegroundColor', col)
end