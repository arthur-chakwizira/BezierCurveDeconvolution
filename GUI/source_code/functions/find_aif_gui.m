function find_aif_gui(handles)
%       This function accepts the handles structure and returns a mean AIF
%       and the area under that mean AIF. It displays a DSC-MRI
%       concentration image on the GUI window and lets the user select
%       individual pixels in search for an AIF.
%        Author:
%              Arthur Chakwizira
%              arthur.chakwizira@med.lu.se
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Note: that handles input is a local copy in this function. If some
%subfunction in here changes handles using guidata(hObject, handles), the
%master copy outside of this function is updated, but not the local copy
%handles. That means subsequent calls to other subfunctions using the local
%handles will be oblivious to the changes performed by other subfunctions.
%To fetch the updated master copy, each subfunction needs to do handles =
%guidata(hObject)
hObject = handles.hObject;

prepare_concentration_data(hObject)
report('AIF Selection: Add ROI to get started', handles)
prepare_sliders(hObject);
set(handles.slider1, 'Callback', {@change_slice});
set(handles.slider8, 'Callback', {@change_timepoint});
initialise_roi_placement(hObject, handles);
initialise_dynamic_aif_plot(hObject, handles);
handles = guidata(hObject);
handles.path_to_aif = [];
%continuously track mouse position on axes1
set(handles.figure1, 'WindowButtonMotionFcn', {@mouseMove});
% set (handles.figure1, 'WindowButtonMotionFcrn', []);
set(handles.pushbutton24, 'Callback', {@clearROIs});
set(handles.pushbutton25, 'Callback', {@saveAIF}, 'String', 'Save AIF');
set(handles.pushbutton30, 'Callback', {@Exit}, 'String', 'Finish AIF Selection');
set(handles.figure1, 'ButtonDownFcn', {@mouseClick})
set(handles.axes1, 'PickableParts', 'none') %to prevent axes and its children from absorbing mouse clicks
set([handles.pushbutton24, handles.pushbutton25, handles.pushbutton30], 'Visible', 'on')
set([handles.pushbutton26], 'Visible', 'off')

%NB! Clear all fields from handles that are not needed elsewhere

handles.data_is_for_aif = false;
guidata(hObject, handles)

end





%% Helper functions
function mouseMove(hObject, ~, ~)
handles = guidata(hObject);
if ~ishandle(handles.aif_conc_image); return; end
C = get(handles.axes1, 'CurrentPoint');
x = round(C(1,1)); %these are flipped
y = round(C(1,2));
sl = get(handles.slider1, 'value');
try
    aif = squeeze(handles.dsc_data_c(y,x,sl,:));
    set(handles.dynamic_aif_plot, 'YData', aif); drawnow
catch
end
end


function prepare_sliders(hObject)
handles = guidata(hObject);
clim_im = [0 100];
aif_slice = ceil(handles.img_size(3)/2);
aif_time_point = handles.t_min_signal;
handles.aif_conc_image = imshow(zeros(handles.img_size(1), handles.img_size(2)), clim_im, 'colormap', gray, 'Parent', handles.axes1);
set(handles.aif_conc_image, 'CData', (squeeze(handles.dsc_data_c(:,:,aif_slice,aif_time_point))))
% title(handles.axes1, ['Concentration: Slice ' num2str(aif_slice) ' : Time point ' num2str(aif_time_point)])
set(handles.slider1, 'Visible', 'on')%vertical scrolling on axes1
set(handles.slider1, 'Min', 1, 'Max', handles.img_size(3), ...
    'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', aif_slice)
set(handles.slider8, 'Visible', 'on')%horizontal scrolling on axes1
set(handles.slider8, 'Min', 1, 'Max', handles.img_size(4), ...
    'SliderStep', [1, 1]/(handles.img_size(4) - 1), 'Value', aif_time_point)
clim_max = double(max(handles.dsc_data_c(isfinite(handles.dsc_data_c)))) + 1E-10;
if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/990; end
if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end %$$$$$$$
set(handles.slider6, 'Visible', 'on') %color setting
set(handles.slider6, 'Min', 10E-10, 'Max', clim_max, ...
    'SliderStep', [step_small, step_big], 'Value', clim_im(2))
set(handles.slider10, 'Visible', 'on') %color setting
set(handles.slider10, 'Min', 0, 'Max', clim_max-1E-10, ...
    'SliderStep', [step_small, step_big], 'Value', 0)
guidata(hObject, handles)
end

function change_slice(hObject, ~, ~)
handles = guidata(hObject);
index = round(get(hObject, 'Value'));
aif_time_point = round(get(handles.slider8, 'Value'));
set(handles.aif_conc_image, 'CData', squeeze(handles.dsc_data_c(:,:,index,aif_time_point)))
if isfield(handles, 'aif_pixel_plots')
    active_rois = handles.aif_pixel_plots(ishandle(handles.aif_pixel_plots));
    active_roi_slices = handles.roi_slices(ishandle(handles.aif_pixel_plots));
    if all(ishandle(handles.aif_pixel_plots)); set(handles.aif_pixel_plots, 'Parent', []);end
    set(active_rois, 'Parent', []);
    set(active_rois(active_roi_slices == index), 'Parent', handles.axes1);
    % title(handles.axes1, ['Concentration: Slice ' num2str(index) ' : Time point ' num2str(aif_time_point)])
end
if isfield(handles, 'aif_pixel_plots'); set(active_rois(active_roi_slices == index), 'Parent', handles.axes1); end
guidata(hObject, handles)
end

function change_timepoint(hObject, ~ , ~)
handles = guidata(hObject);
time_point = round(get(hObject,'Value'));
current_slice = round(get(handles.slider1,'Value'));
set(handles.aif_conc_image, 'CData', squeeze(handles.dsc_data_c(:,:,current_slice, time_point))  )
% title(handles.axes1, ['Concentration' ': Slice ' num2str(current_slice) ' : Time-point ' num2str(time_point)],...
%     'Interpreter', 'None')
guidata(hObject, handles)
end

function initialise_roi_placement(hObject, ~)
handles = guidata(hObject);
max_N_roi = 128*128;
handles.roi_counter = 0;
handles.roi_slices = zeros(max_N_roi, 1);
handles.aif_candidates = NaN(max_N_roi, handles.img_size(4));
handles.opt_aif_coordinates = NaN(max_N_roi, 2);
handles.aif_plots = gobjects(max_N_roi, 1);
handles.aif_pixel_plots = gobjects(max_N_roi, 1);
hold(handles.axes1, 'on')
hold(handles.axes2, 'on')
set(handles.axes1, 'PickableParts', 'none') %to prevent axes and its children from absorbing mouse clicks
handles.t = 0:handles.tr:(handles.img_size(4)-1)*handles.tr;
handles.best_aif_plot = plot(handles.axes2, handles.t, NaN*handles.t, 'k-', 'Linewidth', 2);
guidata(hObject, handles)
end

function initialise_dynamic_aif_plot(hObject, ~)
handles = guidata(hObject);
handles.t = 0:handles.tr:(handles.img_size(4)-1)*handles.tr;
handles.dynamic_aif_plot = plot(handles.axes2,handles.t, handles.t*NaN, 'r:', 'LineWidth', 2);
% legend(handles.axes2, "Mouse tracking")
set(handles.axes2, 'YLim', [-20 120], 'YTick',-20:20:120 , 'XLim', [floor(min(handles.t)) ceil(max(handles.t))],...
    'XTick', floor(min(handles.t)):10:ceil(max(handles.t)), 'box', 'on')
grid(handles.axes2, 'minor')
set(handles.axes2, 'Linewidth', 1.2)
xlabel(handles.axes2, 'time [ms]')
ylabel(handles.axes2, 'Concentration')
legend(handles.axes2, ["Mean AIF", "Mouse tracking"])
guidata(hObject, handles)
end



function prepare_concentration_data(hObject)
handles = guidata(hObject);
if ~isfield(handles, 'dsc_data_c')
    dsc_data_c = handles.dsc_data.*0;
    dsc_data = handles.dsc_data;
    img_size = handles.img_size;
    mask = handles.mask;
    te = handles.te;
    baseline_index = handles.baseline_index;
    for ro_w = 1:img_size(1)
        for co_l = 1:img_size(2)
            for sl_c = 1:img_size(3)
                if mask(ro_w, co_l, sl_c)
                    pixel_signal = squeeze(dsc_data(ro_w, co_l, sl_c, :));
                    pixel_s0 = mean(pixel_signal(baseline_index));
                    pixel_conc = (-1/te)*log(pixel_signal./pixel_s0);
                    dsc_data_c(ro_w, co_l, sl_c,:) = pixel_conc;
                else
                    dsc_data_c(ro_w, co_l, sl_c,:) = 0;
                end
            end
        end
    end
    dsc_data_c(~isfinite(dsc_data_c)) = 0;
    dsc_data_c(isnan(dsc_data_c)) = 0;
    handles.dsc_data_c = dsc_data_c;
    clear('dsc_data_c', 'dsc_data', 'mask')
    guidata(hObject, handles)
end
if ishandle(handles.wait_for_data); delete(handles.wait_for_data); end
end

function clearROIs(hObject, ~, ~)
%clear all drawn ROIS
handles = guidata(hObject);
% set(handles.axes1, 'PickableParts', 'all') %to prevent axes and its children from absorbing mouse clicks
clim_im = get(handles.axes1, 'CLim');
active_rois = handles.aif_pixel_plots(ishandle(handles.aif_pixel_plots));
% active_roi_slices = handles.roi_slices(ishandle(handles.rois));
for c_ar = 1:numel(active_rois)
    delete(active_rois(c_ar))
end
if isfield(handles, 'aif_plots'); delete(handles.aif_plots); handles = rmfield(handles, 'aif_plots'); end
if isfield(handles, 'best_aif_plot'); delete(handles.best_aif_plot); handles = rmfield(handles, 'best_aif_plot'); end
if isfield(handles, 'aif_pixel_plots'); delete(handles.aif_pixel_plots); handles = rmfield(handles, 'aif_pixel_plots'); end
current_slice = round(get(handles.slider1, 'Value'));
current_timepoint = round(get(handles.slider8, 'Value'));
cla(handles.axes1, 'reset')
handles.aif_conc_image = imshow(zeros(handles.img_size(1), handles.img_size(2)), clim_im, 'colormap', gray, 'Parent', handles.axes1);
set(handles.aif_conc_image, 'CData', (squeeze(handles.dsc_data_c(:,:,current_slice,current_timepoint))) )
cla(handles.axes2, 'reset')
guidata(hObject, handles)
initialise_roi_placement(hObject)
initialise_dynamic_aif_plot(hObject)
% handles = guidata(hObject);
end

function saveAIF(hObject, ~, ~)
handles = guidata(hObject);
[file, path] = uiputfile( '*.*','Save AIF', fullfile(handles.file_folder, strcat('aif_',handles.file_name)));
if file == 0
    report('AIF not saved.', handles, 'r')
else
    path_to_aif = fullfile(path,file);
    handles.path_to_aif = path_to_aif;
    guidata(hObject, handles)
    save_data.header_info = [];
    save_data.this_folder = path_to_aif; save_data.this_format = 0;
    save_data.this_name = 0; save_data.data_to_save = handles.best_aif; save_this_file(save_data)
end
end

function mouseClick(hObject, ~, ~)
handles = guidata(hObject);
% if ~ishandle(handles.aif_conc_image); return; end
% disp("Clicked")
hold(handles.axes1, 'on')
C = get(handles.axes1, 'CurrentPoint');
x = round(C(1,1)); %these are flipped
y = round(C(1,2));
pixel_x = x;
pixel_y = y;

if ~ismember(pixel_x, 1:handles.img_size(2)); return; end
if ~ismember(pixel_y, 1:handles.img_size(1)); return; end


current_slice = round(get(handles.slider1, 'Value'));
this_indx = false;
check_exists = (handles.roi_slices == current_slice & handles.opt_aif_coordinates(:,1) == pixel_x & ...
    handles.opt_aif_coordinates(:, 2) == pixel_y);
exist_s = any(check_exists);
if exist_s; this_indx = find(check_exists); end

if exist_s && this_indx ~= false %check if the clicked pixel is one of our optimal pixels
    data_on_pixel = get(handles.aif_pixel_plots(this_indx), 'YData'); %check if its plot exists
    if ~isempty(data_on_pixel) %if it does
        set(handles.aif_pixel_plots(this_indx), 'XData', [], 'YData', []) %delete it
        set( handles.aif_plots(this_indx), 'XData', [], 'YData', []) %and delete the corresponding aif curve
        handles.aif_candidates(this_indx, :) = NaN(1, handles.img_size(4)); %remove the curve from the list of optimal ones
        handles.best_aif = nanmean(handles.aif_candidates, 1); %update the mean
        set(handles.best_aif_plot, 'XData', handles.t, 'YData', handles.best_aif) %and update the mean plot
        drawnow     %now
    else %if the plot nolonger exists (already clicked off)
        set(handles.aif_pixel_plots(this_indx), 'XData', handles.opt_aif_coordinates(this_indx, 1), 'YData', handles.opt_aif_coordinates(this_indx, 2)) %restore it
        set( handles.aif_plots(this_indx), 'XData', handles.t, 'YData', handles.dsc_data_c(pixel_y, pixel_x, current_slice, :)) %restore its aif curve
        handles.aif_candidates(this_indx, :) = handles.dsc_data_c(pixel_y, pixel_x, current_slice, :); %restore the curve to the list of optimal ones
        handles.best_aif = nanmean(handles.aif_candidates, 1); %update the mean
        set(handles.best_aif_plot, 'XData', handles.t, 'YData', handles.best_aif) %and update the mean plot
        drawnow     %now
    end
else %if clicked pixel is not one of our optimal ones
    handles.roi_counter = handles.roi_counter +1;
    handles.aif_candidates(handles.roi_counter, :) = handles.dsc_data_c(pixel_y, pixel_x, current_slice, :);
    handles.best_aif = nanmean(handles.aif_candidates, 1); %update the mean
    handles.opt_aif_coordinates(handles.roi_counter, :) = [pixel_x, pixel_y];
    handles.roi_slices(handles.roi_counter) = current_slice;
    set(handles.best_aif_plot, 'XData', handles.t, 'YData', handles.best_aif) %and update the mean plot
    handles.aif_plots(handles.roi_counter) = plot(handles.axes2, handles.t,handles.aif_candidates(handles.roi_counter, :), 'Linewidth', 1, 'HandleVisibility', 'off');
    col = handles.aif_plots(handles.roi_counter).Color;
    handles.aif_pixel_plots(handles.roi_counter) = plot(handles.axes1, pixel_x, pixel_y, 's', 'MarkerFaceColor', col, 'MarkerEdgeColor', col,...
        'MarkerSize', 5);
    drawnow
end % add it to the collection
handles.aif_area = trapz(handles.t, handles.best_aif);
guidata(hObject, handles)
end

function Exit(hObject, ~, ~)
handles = guidata(hObject);
try
    handles = rmfield(handles, {'aif_candidates', 'dsc_data_c', 'aif_plots', 'aif_pixel_plots', 'opt_aif_coordinates', 'roi_slices'});
catch
end
% set(handles.slider1, 'Callback', handles.slider1_Callback);
% set(handles.slider8, 'Callback', handles.slider8_Callback);
%continuously track mouse position on axes1
set(handles.slider1, 'Callback', []);
set(handles.slider8, 'Callback', []);
set(handles.figure1, 'WindowButtonMotionFcn', []);
% set (handles.figure1, 'WindowButtonMotionFcrn', []);
set(handles.pushbutton24, 'Callback', []);
set(handles.pushbutton25, 'Callback', []);
set(handles.pushbutton30, 'Callback',[]);
set(handles.figure1, 'ButtonDownFcn',[])
set(handles.axes1, 'PickableParts', 'all') %to prevent axes and its children from absorbing mouse clicks
reset_axes(handles)
report('AIF selection closed.', handles)
set([handles.pushbutton24, handles.pushbutton25, handles.pushbutton30], 'Visible', 'off')
guidata(hObject, handles)
end