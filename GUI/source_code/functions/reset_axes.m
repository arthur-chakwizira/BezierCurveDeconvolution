function reset_axes(handles, ax)
%restore axes to original appearance when GUI was opened
if nargin < 2
cla(handles.axes1, 'reset')
cla(handles.axes2, 'reset')
set(handles.slider1, 'Visible', 'off') %slice change on axes1
set(handles.slider2, 'Visible', 'off'); %slice change on axes2
set(handles.slider6, 'Visible', 'off') %upper window setting on axes1
set(handles.slider10, 'Visible', 'off')%lower window setting on axes1
set(handles.slider7, 'Visible', 'off')%upper window setting on axes2
set(handles.slider11, 'Visible', 'off')%lower window setting on axes2
set(handles.slider8, 'Visible', 'off')%horizontal scrolling on axes1
set(handles.slider9, 'Visible', 'off')%horizontal scrolling on axes2
set(handles.axes1,'XTick', [], 'YTick', []) %set default appearance for axes
set(handles.axes2,'XTick', [], 'YTick', [])
set(handles.axes1, 'box', 'on')
set(handles.axes2, 'box', 'on')
set([handles.pushbutton15, handles.pushbutton16, handles.pushbutton17, handles.pushbutton18, ...
    handles.togglebutton2, handles.togglebutton3], 'Visible', 'off')
else
cla(ax, 'reset')
set(ax,'XTick', [], 'YTick', []) %set default appearance for axes
set(ax, 'box', 'on')
if ax == handles.axes1
set([handles.pushbutton15, handles.pushbutton17, ...
    handles.togglebutton2], 'Visible', 'off') 
set([handles.slider1, handles.slider8, handles.slider6, handles.slider10], 'Visible', 'off')
end

if ax == handles.axes2
set([handles.pushbutton16, handles.pushbutton18, ...
    handles.togglebutton3], 'Visible', 'off') 
set([handles.slider2, handles.slider9, handles.slider7, handles.slider11], 'Visible', 'off')
end

end

set(handles.figure1, 'WindowButtonMotionFcn', []);
end