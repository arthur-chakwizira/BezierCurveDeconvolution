function terminate_all_aif_or_vof_selection(hObject)
%This needs to be called on "Load DSC data"
Exit(hObject)
% disp("Successfully stopped AIF selection")
end

function Exit(hObject, ~, ~)
handles = guidata(hObject);
try
handles = rmfield(handles, {'aif_candidates', 'aif_plots', 'aif_pixel_plots', 'opt_aif_coordinates', 'roi_slices', 'rois', 'best_aif'});
catch
%     disp("Some fields could not be removed from handles")
end
try handles = rmfield(handles,'dsc_data_c'); catch; end
set(handles.slider1, 'Callback', []);
set(handles.slider8, 'Callback', []);
%continuously track mouse position on axes1
set(handles.figure1, 'WindowButtonMotionFcn', []);
% set (handles.figure1, 'WindowButtonMotionFcrn', []);
set(handles.pushbutton24, 'Callback', []);
set(handles.pushbutton25, 'Callback', []);
set(handles.pushbutton26, 'Callback', []);
set(handles.pushbutton27, 'Callback', []);
set(handles.pushbutton30, 'Callback',[]);
set(handles.figure1, 'ButtonDownFcn',[])
set(handles.axes1, 'PickableParts', 'all') %to prevent axes and its children from absorbing mouse clicks
reset_axes(handles)
report('AIF selection closed.', handles)
set([handles.pushbutton24, handles.pushbutton25,handles.pushbutton26, handles.pushbutton27, handles.pushbutton30], 'Visible', 'off')
guidata(hObject, handles)
end
