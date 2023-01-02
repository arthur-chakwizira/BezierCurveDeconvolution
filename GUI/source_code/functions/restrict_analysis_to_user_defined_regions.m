function restrict_analysis_to_user_defined_regions(hObject)
%allows user to select some ROIs in which to perform analysis
handles = guidata(hObject);

set([handles.pushbutton31, handles.pushbutton32,handles.pushbutton33, handles.pushbutton34, handles.pushbutton35], 'Visible', 'off') %initially turn off all these buttons
set([handles.pushbutton24, handles.pushbutton25,handles.pushbutton26, handles.pushbutton27, handles.pushbutton30], 'Visible', 'off') %initially turn off all these buttons
reset_axes(handles, handles.axes1)

set(handles.pushbutton30, 'String', 'Continue');
set(handles.pushbutton30, 'Visible', 'on');
set(handles.pushbutton30, 'Callback', {@Proceed})
set(handles.pushbutton26, 'Visible', 'on');
set(handles.pushbutton24, 'Visible', 'on');
set(handles.pushbutton26, 'Callback', {@addROI});
set(handles.pushbutton24, 'Callback', {@clearROIs});
set(handles.slider1, 'Callback', {@change_slice});
set(handles.slider8, 'Callback', {@change_timepoint});
prepare_dsc_data(hObject)
prepare_sliders(hObject)
initialise_roi_placement(hObject)
waitfor(handles.pushbutton30, 'Visible')
handles = guidata(hObject);
active_rois = handles.rois(ishandle(handles.rois));
active_roi_slices = handles.roi_slices(ishandle(handles.rois));
new_mask = false(handles.img_size(1), handles.img_size(2),  handles.img_size(3));

img_size = handles.img_size;
for c_roi = 1:numel(active_rois) %use meshgrid if this is slow
    h_roi = active_rois(c_roi);
    oef_slice = active_roi_slices(c_roi);
    for x = 1:img_size(1)
        for y = 1:img_size(2)
            if inROI(h_roi, x,y)
                new_mask(y,x,oef_slice) = true;
            end
        end
    end
end

%update mask in handles
%new mask is true only where old mask was also true
handles.old_mask = handles.mask; %keep this for later
handles.mask = handles.mask & new_mask;
guidata(hObject, handles)

handles = rmfield(handles, 'dsc_data_filtered'); %remove this to save memory
guidata(hObject, handles)
end


%% Helper functions
function prepare_sliders(hObject, ~, ~)
handles = guidata(hObject);
initial_slice = round(handles.img_size(3)/2);
t_min_signal = handles.t_min_signal;

initial_slice_image = squeeze(handles.dsc_data_filtered(:,:,initial_slice,handles.t_min_signal));
clim_im = [0   max(initial_slice_image(isfinite(initial_slice_image)))]; %color limit, can be changed using sliders
handles.dsc_data_image = imshow(zeros(handles.img_size(1), handles.img_size(2)),clim_im, 'colormap', gray, 'Parent', handles.axes1);
% title(handles.axes1, ['DSC signal: Slice ' num2str(initial_slice), ' : Time point ' num2str(handles.t_min_signal)])
set(handles.dsc_data_image, 'CData', initial_slice_image ) %display initial slice
set(handles.slider1, 'Min', 1, 'Max', handles.img_size(3), ... %setup slider for changing slice
    'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', initial_slice)
set(handles.slider8, 'Visible', 'on')%horizontal scrolling on axes1
set(handles.slider1, 'Visible', 'on')%horizontal scrolling on axes1
set(handles.slider8, 'Min', 1, 'Max', handles.img_size(4), ...
    'SliderStep', [1, 1]/(handles.img_size(4) - 1), 'Value', t_min_signal)

clim_max = double(max(handles.dsc_data_filtered(isfinite(handles.dsc_data_filtered)))) + 1E-10;
if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/999; end
if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end %setup slider for changing window setting
set(handles.slider6, 'Visible', 'on') %color setting
set(handles.slider6, 'Min', 10e-10, 'Max', clim_max, ...
    'SliderStep', [step_small, step_big], 'Value', clim_im(2))
set(handles.slider10, 'Visible', 'on') %color setting
set(handles.slider10, 'Min', 0, 'Max', clim_max-1E-10, ...
    'SliderStep', [step_small, step_big], 'Value', 0)
% hObject = handles.output;
title(handles.axes1, ['Slice ' num2str(initial_slice)])
guidata(hObject, handles)
end

function initialise_roi_placement(hObject, ~)
handles = guidata(hObject);
max_N_roi = 128*128;
handles.roi_counter = 0;
handles.rois = gobjects(max_N_roi,1);
handles.roi_slices = zeros(max_N_roi, 1);
handles.roi_mask = false(handles.img_size(1), handles.img_size(2), handles.img_size(3));
guidata(hObject, handles)
end

function Proceed(hObject, ~, ~)
%this is called when ROI placement is complete
handles = guidata(hObject);
set(handles.pushbutton26, 'Visible', 'off');
set(handles.pushbutton24, 'Visible', 'off');
set(handles.pushbutton24, 'Callback', [])
set(handles.pushbutton26, 'Callback', [])
set(handles.pushbutton30, 'String', 'Finish AIF selection')
set(handles.pushbutton30, 'Callback', [])
set(handles.slider1, 'Callback', []);
set(handles.slider8, 'Callback', []);
set(handles.pushbutton30, 'Visible', 'off')
guidata(hObject, handles)
end


function addROI(hObject, ~, ~)
%adds 1 ROI when called
handles = guidata(hObject);
roi = drawfreehand(handles.axes1, 'Color', 'r', 'LineWidth',1, 'FaceAlpha', 0);
handles.roi_counter = handles.roi_counter+1;
handles.rois(handles.roi_counter) = roi;
handles.roi_slices(handles.roi_counter) = get(handles.slider1, 'value');
guidata(hObject, handles) %update master copy
end


function clearROIs(hObject, ~, ~)
%clear all drawn ROIS
handles = guidata(hObject);
clim_im = get(handles.axes1, 'CLim');
active_rois = handles.rois(ishandle(handles.rois));
% active_roi_slices = handles.roi_slices(ishandle(handles.rois));
for c_ar = 1:numel(active_rois)
    delete(active_rois(c_ar))
end
current_slice = round(get(handles.slider1, 'Value'));
current_timepoint = round(get(handles.slider8, 'Value'));
cla(handles.axes1, 'reset')
handles.dsc_data_image = imshow(zeros(handles.img_size(1), handles.img_size(2)), clim_im, 'colormap', gray, 'Parent', handles.axes1);
set(handles.dsc_data_image, 'CData', (squeeze(handles.dsc_data_filtered(:,:,current_slice,current_timepoint))) )
guidata(hObject, handles)
initialise_roi_placement(hObject)
end


function change_slice(hObject, ~, ~)
handles = guidata(hObject);
current_slice = round(get(hObject, 'Value'));
current_time_point = round(get(handles.slider8, 'Value'));
set(handles.dsc_data_image, 'CData', squeeze(handles.dsc_data_filtered(:,:,current_slice,current_time_point)))
active_rois = handles.rois(ishandle(handles.rois));
active_roi_slices = handles.roi_slices(ishandle(handles.rois));
set(active_rois, 'Parent', []);
set(active_rois(active_roi_slices == current_slice), 'Parent', handles.axes1);
title(handles.axes1, ['Slice ' num2str(current_slice)])
guidata(hObject, handles)
end


function change_timepoint(hObject, ~ , ~)
handles = guidata(hObject);
time_point = round(get(hObject,'Value'));
current_slice = round(get(handles.slider1,'Value'));
set(handles.dsc_data_image, 'CData', squeeze(handles.dsc_data_filtered(:,:,current_slice, time_point))  )
title(handles.axes1, ['Slice ' num2str(current_slice)])
% title(handles.axes1, ['Signal' ': Slice ' num2str(current_slice) ' : Time-point ' num2str(time_point)],...
%     'Interpreter', 'None')
guidata(hObject, handles)
end

function prepare_dsc_data(hObject)
handles = guidata(hObject);
if isfield(handles, 'mask') %if mask has been generated
    handles.dsc_data_filtered = handles.dsc_data.*0;
    for sl = 1:handles.img_size(3) %filter image to be displayed using the mask
        for tp = 1:handles.img_size(4)
            handles.dsc_data_filtered(:,:,sl, tp) = squeeze(handles.dsc_data(:,:,sl, tp)).*squeeze(handles.mask(:,:,sl));
        end
    end
else
    handles.dsc_data_filtered = handles.dsc_data;
end
guidata(hObject, handles)
end