function [mean_vof_c, vof_area] = find_vof_semi_auto_gui(handles)
%       This function accepts the handles structure and returns a mean VOF
%       and the area under that mean VOF. It displays a DSC-MRI
%       concentration image on the GUI window and lets the user place a
%       ROI. The function suggests VOFs in the drawn ROI and allows the
%       user to edit the suggestion.
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
report('VOF Selection: Add ROI to get started', handles)
prepare_sliders(hObject);
set(handles.slider1, 'Callback', {@change_slice});
set(handles.slider8, 'Callback', {@change_timepoint});
initialise_roi_placement(hObject, handles);
initialise_dynamic_vof_plot(hObject, handles);
handles = guidata(hObject);
handles.path_to_vof = [];
%continuously track mouse position on axes1
set(handles.figure1, 'WindowButtonMotionFcn', {@mouseMove});
% set (handles.figure1, 'WindowButtonMotionFcrn', []);
set(handles.pushbutton26, 'Callback', {@addROI});
set(handles.pushbutton27, 'Callback', {@suggest_vofs}, 'String', 'Suggest VOFs');
set(handles.pushbutton24, 'Callback', {@clearROIs});
set(handles.pushbutton25, 'Callback', {@saveVOF}, 'String', 'Save VOF');
set(handles.pushbutton30, 'Callback', {@Exit}, 'String', 'Finish VOF Selection');
set([handles.pushbutton24, handles.pushbutton25,handles.pushbutton26, handles.pushbutton27, handles.pushbutton30], 'Visible', 'on')

set(handles.figure1, 'ButtonDownFcn', {@mouseClick})
mean_vof_c = 0;
vof_area = 0;

%NB! Clear all fields from handles that are not needed elsewhere

handles.data_is_for_vof = false;
guidata(hObject, handles)

end





%% Helper functions
function mouseMove(hObject, ~, ~)
handles = guidata(hObject);
if ~ishandle(handles.vof_conc_image); return; end
C = get(handles.axes1, 'CurrentPoint');
x = round(C(1,1)); %these are flipped
y = round(C(1,2));
sl = get(handles.slider1, 'value');
try
    vof = squeeze(handles.dsc_data_c(y,x,sl,:));
    set(handles.dynamic_vof_plot, 'YData', vof); drawnow
catch
end
end


function suggest_vofs(hObject, ~, ~)
% tic
handles = guidata(hObject);
if ~ishandle(handles.vof_conc_image)
    report('There was no image on axes1. Draw ROI to continue', handles)
    cla(handles.axes1, 'reset')
    prepare_sliders(hObject)
    clearROIs(hObject)
    return
end

img_size = handles.img_size;
mask = handles.mask;
dsc_data_c = handles.dsc_data_c;
dsc_data = handles.dsc_data;
vof_mask = false(img_size(1:3));
baseline_idx = handles.baseline_index;
tail_idx = handles.tail_index;
t = handles.t;

cla(handles.axes2, 'reset');
initialise_dynamic_vof_plot(hObject);
handles = guidata(hObject);

if isfield(handles, 'vof_plots')
for c_p = 1:numel(handles.vof_plots)
    if ishandle(handles.vof_plots(c_p)); delete(handles.vof_plots(c_p));end
    if ishandle(handles.vof_pixel_plots(c_p)); delete(handles.vof_pixel_plots(c_p)); end
end
if ishandle(handles.best_vof_plot); delete(handles.best_vof_plot);end
end

% %go through all active drawn ROIs
active_rois = handles.rois(ishandle(handles.rois));
active_roi_slices = handles.roi_slices(ishandle(handles.rois));
hold(handles.axes1, 'on')

%initialise candidate vof plots
hold(handles.axes2, 'on')
handles.vof_plots = gobjects(5,1);
handles.vof_pixel_plots = gobjects(5,1);
handles.vof_plot_slices = zeros(5,1);
for c_vof_plot = 1:numel(handles.vof_plots)
    handles.vof_plots(c_vof_plot) = plot(handles.axes2, t, t*NaN, '-', 'Linewidth', 1,...
        'Handlevisibility', 'off', 'Parent', handles.axes2);
    handles.vof_pixel_plots(c_vof_plot) = plot(handles.axes2, NaN, NaN, '+', 'Linewidth', 1,...
        'Handlevisibility', 'off', 'Parent',  handles.axes1);
end
%initialise best vof plot
handles.best_vof_plot = plot(handles.axes2, t, t*NaN, 'k-', 'Linewidth', 2);
drawnow

%extract coordinates of all pixels in the drawn ROIs
pixel_xcoords = NaN(img_size(1)*img_size(2)*numel(active_rois),1);
pixel_ycoords = NaN(img_size(1)*img_size(2)*numel(active_rois),1);
pixel_slices = NaN(img_size(1)*img_size(2)*numel(active_rois),1);

this_pos = 0;
for c_roi = 1:numel(active_rois) %use meshgrid if this is slow
    h_roi = active_rois(c_roi);
    vof_slice = active_roi_slices(c_roi);
    for x = 1:img_size(1)
        for y = 1:img_size(2)
            if inROI(h_roi, x,y)
                this_pos = this_pos + 1;
                pixel_xcoords(this_pos) = x;
                pixel_ycoords(this_pos) = y;
                pixel_slices(this_pos) = vof_slice;
            end
        end
    end
end
%delete nan values
pixel_xcoords((this_pos+1):end) = [];
pixel_ycoords((this_pos+1):end) = [];
pixel_slices((this_pos+1):end) = [];

%forbid point ROIs
if (isempty(pixel_xcoords) || numel(pixel_xcoords) < 5) %@@@@@
    report('Drawn ROIs contain no useful pixels. Try larger ROI.', handles)
else
    
    vof_c = NaN(size(pixel_xcoords, 1),  img_size(4));
    vof_s = NaN(size(pixel_xcoords, 1),  img_size(4));
    
    for k = 1:length(pixel_xcoords)
        vof_slice = pixel_slices(k);
        vof_s(k,:) = squeeze(dsc_data(pixel_ycoords(k), pixel_xcoords(k), vof_slice, :));
        vof_c(k,:) = squeeze(dsc_data_c(pixel_ycoords(k), pixel_xcoords(k), vof_slice, :));
    end
    
    disqualified = false(numel(pixel_xcoords), 1);
    for ii = 1:numel(pixel_xcoords) %loop through all pixels in ROI
        [vof_max, tmax] = max(vof_c(ii,:)); %maximum of VOF concentration and time point at which it occurs
        vof_baseline = mean(vof_s(ii, baseline_idx)); %baseline of VOF signal
        vof_tail = mean(vof_s(ii, tail_idx)); %tail of the VOF signal
        t20percent = find(vof_c(ii,:) > 0.2*vof_max); %find where VOF concentration is above 20% of its maximum (that is, locate peak)
        if isempty(t20percent) %happens if curve is just a bunch of zeros
            disqualified(ii) = true;  %discriminate against this
        else
            num_peaks_check = diff(t20percent); %checking for curves with multiple distinct maxima
            if ~all(num_peaks_check == 1) %if there is a jump of more than 1 in the peak indices, the curve has more than one peak
                disqualified(ii) = true;
            else %if there is just one peak (it may still have small subpeaks that are not detected by the check above)
                t20_low = t20percent(1); %extract beginning
                t20_up = t20percent(end); %and end of peak
                rising_end = diff(vof_c(ii, t20_low:tmax)); %check slope on rising side
                falling_end = diff(vof_c(ii, tmax:t20_up));  %check slope on falling side
                if ~all(rising_end > 0); disqualified(ii) = true; end %if not monotonically increasing, discriminate
                if ~all(falling_end < 0); disqualified(ii) = true; end %if not monotonically decreasing, discriminate
            end
        end
        %     if vof_baseline < 0.6*max_baseline; vof_c(ii, :) = 0; end %if VOF signal baseline is lower than 60% of maximum in ROI, discriminate
        %     if vof_tail < 0.8*vof_baseline; vof_c(ii, :) = 0; end %if VOF signal tail is lower than 80% of original baseline, discriminate
        %     if min(vof_s(ii,:)) > 0.2*vof_baseline; disqualified(ii) = true; end %if VOF signal drop is lower than 80% of baseline, discriminate
        if any(vof_c(ii, :) < -10); disqualified(ii) = true; end
        if any(isnan(vof_c(ii,:)));  disqualified(ii) = true; end
    end
    
    pixel_xcoords(disqualified) = [];
    pixel_ycoords(disqualified) = [];
    pixel_slices(disqualified) = [];
    vof_c(disqualified, :) = [];
    vof_s(disqualified, :) = [];
    
    
    ymax    = zeros(length(pixel_xcoords), 1);
    fm      = zeros(length(pixel_xcoords), 1);
    noise   = zeros(length(pixel_xcoords), 1);
    pointy  = zeros(length(pixel_xcoords), 1);
    
    vof_data_c = vof_c;
    for i = 1:length(pixel_xcoords)
        [ymax(i), fm(i)] = max(vof_data_c(i,:)); %extract amplitude and time of occurence
        %     fm(i) = sum(t.*vof_data_c(i,:))/sum(vof_data_c(i,:));
        %fm(i) = sum(t(fit_idx).*vof_data_c(i,fit_idx))/sum(vof_data_c(i,fit_idx));
        noise(i) = std(vof_data_c(i,baseline_idx)) + std(vof_data_c(i,tail_idx)); %compute noise level (baseline + tail)
        pointy(i) = sum(abs(diff(vof_data_c(i,:)))); %compute pointyness
    end
    
    numvofs = 5; %number of curves to work with
    if numvofs > length(pixel_xcoords); numvofs = length(pixel_xcoords); end %if greater than number of selected pixels, set to number of selected pixels
    [~,sortIndex] = sort(ymax(:),'descend'); %sort in descending order
    ymax_index = sortIndex(1:numvofs); % get highest 10
    
    % [~,sortIndex] = sort(noise(:),'ascend'); %sort in ascending noise
    % noise_index = sortIndex(1:numvofs); % get lowest 10
    %
    % [~,sortIndex] = sort(pointy(:),'ascend'); % sort in scending pointyness
    % pointy_index = sortIndex(1:numvofs); %get lowest 10
    %
    [~,sortIndex] = sort(fm(:),'descend'); %sort in ascending fm
    fm_index = sortIndex(1:numvofs); % get lowest 10
    
    
    % % find intersects
    % ymaxfm_index = intersect(intersect(intersect(ymax_index, fm_index), noise_index), pointy_index); %initially demand satisfaction of all requirements
    % if (isempty(ymaxfm_index)||length(ymaxfm_index)<5); ymaxfm_index =intersect(intersect(ymax_index, fm_index), noise_index); end %if failed, remove low pointyness requirement
    % if (isempty(ymaxfm_index)||length(ymaxfm_index)<5); ymaxfm_index = intersect(ymax_index, fm_index); end %if still no candidates, remove noise requirement
    ymaxfm_index = intersect(ymax_index, fm_index);
    numvofs = length(ymaxfm_index); %check how many VOFs we have left
    if numvofs < 5 %if none
        ymaxfm_index = ymax_index; %use only those with earliest arrival, that is, remove all other requirements
        numvofs = length(ymaxfm_index); %check how many we have left now
    end
    if numvofs > 5; numvofs = 5; ymaxfm_index = ymaxfm_index(1:5); end %if more than 5, pick only the first 5
    
    %make this change: //17 nov 2022
    % ymaxfm_index = ymax_index;
    
    %plot curves deemed optimal
    hold(handles.axes1, 'on')
    title(handles.axes2, 'VOF selection')
    legend(handles.axes2, ["Mouse-tracking", "Mean VOF"])
    hold(handles.axes2, 'on')
    opt_vof_c = NaN(numvofs, img_size(4));
    handles.opt_vof_coordinates = zeros(numvofs, 2); %x y
    handles.opt_vof_slices = zeros(numvofs, 1);
    %     opt_colors = zeros(numvofs,1);
    
    color_list = {[0 0 0.8095], [1 0.3858 0], [0.8254 0 1], [0.65 0. 0], [0 0.8189 0]};
    for k = 1:numvofs
        opt_vof_c(k,:) = vof_data_c(ymaxfm_index(k), :);
        handles.opt_vof_slices(k) = pixel_slices(ymaxfm_index(k));
        set(handles.vof_plots(k), 'YData', opt_vof_c(k, :), 'Parent', handles.axes2, 'Color', color_list{k});
        this_color = get(handles.vof_plots(k), 'color');
        %         opt_colors(k) = this_color; %we will save these for later
        if all(isnan(opt_vof_c(k, :))); marker_size = eps; else; marker_size = 12; end %@@@
        set(handles.vof_pixel_plots(k), 'XData', pixel_xcoords(ymaxfm_index(k)), 'YData', pixel_ycoords(ymaxfm_index(k)), ...
            'MarkerFaceColor', this_color, 'MarkerEdgeColor', this_color, 'MarkerSize', marker_size, 'LineWidth', 2, 'Parent', []);
        handles.opt_vof_coordinates(k,:) = [pixel_xcoords(ymaxfm_index(k)),   pixel_ycoords(ymaxfm_index(k))];
    end
    
    current_slice = round(get(handles.slider1, 'Value'));
    set(handles.vof_pixel_plots(handles.opt_vof_slices == current_slice), 'Parent', handles.axes1)
    %make this change: //17 nov 2022
    % ymaxfm_index = ymax_index;
    
    candidates = vof_data_c(ymaxfm_index,:);
    best_vof = mean(candidates, 1); %take average.
    set(handles.best_vof_plot, 'YData', best_vof);
    drawnow
    %     uistack(handles.dynamic_vof_plot, 'top')
    handles.vof_area = trapz(t,best_vof);
    handles.best_vof = best_vof;
    handles.vof_candidates = candidates;
    handles.numvofs = numvofs;
    handles.ready_to_edit = true; %used to inform edit function that no clicks have been made yet
end

actice_rois = handles.rois(ishandle(handles.rois)); %make space for editing vof suggestions
for c_ar = 1:numel(actice_rois)
    set(actice_rois(c_ar), 'LineWidth', 0.5, 'Waypoints', false(size(actice_rois(c_ar).Waypoints)), 'FaceAlpha', 0, 'FaceSelectable', 0)
end
set(handles.axes1, 'PickableParts', 'none') %to prevent axes and its children from absorbing mouse clicks
% set(handles.vof_conc_image, 'ButtonDownFcn', {@mouseClick})
guidata(hObject, handles)
% toc
end



function addROI(hObject, ~, ~)
handles = guidata(hObject);
set(handles.axes1, 'PickableParts', 'all') %to prevent axes and its children from absorbing mouse clicks
roi = drawfreehand(handles.axes1, 'Color', 'r', 'LineWidth',1, 'FaceAlpha', 0);
handles.roi_counter = handles.roi_counter+1;
handles.rois(handles.roi_counter) = roi;
handles.roi_slices(handles.roi_counter) = get(handles.slider1, 'value');
guidata(hObject, handles) %update master copy
end


function prepare_sliders(hObject)
handles = guidata(hObject);
clim_im = [0 100];
vof_slice = ceil(handles.img_size(3)/2);
vof_time_point = handles.t_min_signal;
handles.vof_conc_image = imshow(zeros(handles.img_size(1), handles.img_size(2)), clim_im, 'colormap', gray, 'Parent', handles.axes1);
set(handles.vof_conc_image, 'CData', (squeeze(handles.dsc_data_c(:,:,vof_slice,vof_time_point))))
% title(handles.axes1, ['Concentration: Slice ' num2str(vof_slice) ' : Time point ' num2str(vof_time_point)])
set(handles.slider1, 'Visible', 'on')%vertical scrolling on axes1
set(handles.slider1, 'Min', 1, 'Max', handles.img_size(3), ...
    'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', vof_slice)
set(handles.slider8, 'Visible', 'on')%horizontal scrolling on axes1
set(handles.slider8, 'Min', 1, 'Max', handles.img_size(4), ...
    'SliderStep', [1, 1]/(handles.img_size(4) - 1), 'Value', vof_time_point)
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
vof_time_point = round(get(handles.slider8, 'Value'));
set(handles.vof_conc_image, 'CData', squeeze(handles.dsc_data_c(:,:,index,vof_time_point)))
active_rois = handles.rois(ishandle(handles.rois));
active_roi_slices = handles.roi_slices(ishandle(handles.rois));
if isfield(handles, 'vof_pixel_plots'); if all(ishandle(handles.vof_pixel_plots)); set(handles.vof_pixel_plots, 'Parent', []);end; end
set(active_rois, 'Parent', []);
set(active_rois(active_roi_slices == index), 'Parent', handles.axes1);
% title(handles.axes1, ['Concentration: Slice ' num2str(index) ' : Time point ' num2str(vof_time_point)])
if isfield(handles, 'vof_pixel_plots')&&isfield(handles, 'opt_vof_slices'); set(handles.vof_pixel_plots(handles.opt_vof_slices == index), 'Parent', handles.axes1); end
guidata(hObject, handles)
end

function change_timepoint(hObject, ~ , ~)
handles = guidata(hObject);
time_point = round(get(hObject,'Value'));
current_slice = round(get(handles.slider1,'Value'));
set(handles.vof_conc_image, 'CData', squeeze(handles.dsc_data_c(:,:,current_slice, time_point))  )
% title(handles.axes1, ['Concentration' ': Slice ' num2str(current_slice) ' : Time-point ' num2str(time_point)],...
%     'Interpreter', 'None')
end

function initialise_roi_placement(hObject, ~)
handles = guidata(hObject);
max_N_roi = 128*128;
handles.roi_counter = 0;
handles.rois = gobjects(max_N_roi,1);
handles.roi_slices = zeros(max_N_roi, 1);
guidata(hObject, handles)
end

function initialise_dynamic_vof_plot(hObject, ~)
handles = guidata(hObject);
handles.t = 0:handles.tr:(handles.img_size(4)-1)*handles.tr;
handles.dynamic_vof_plot = plot(handles.axes2,handles.t, handles.t*NaN, 'r:', 'LineWidth', 2);
set(handles.axes2, 'YLim', [-20 120], 'YTick',-20:20:120 , 'XLim', [floor(min(handles.t)) ceil(max(handles.t))])
grid(handles.axes2, 'minor')
set(handles.axes2, 'Linewidth', 1.2)
xlabel(handles.axes2, 'time [ms]')
ylabel(handles.axes2, 'Concentration')
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
set(handles.axes1, 'PickableParts', 'all') %to prevent axes and its children from absorbing mouse clicks
clim_im = get(handles.axes1, 'CLim');
active_rois = handles.rois(ishandle(handles.rois));
% active_roi_slices = handles.roi_slices(ishandle(handles.rois));
for c_ar = 1:numel(active_rois)
    delete(active_rois(c_ar))
end
if isfield(handles, 'vof_plots'); delete(handles.vof_plots); handles = rmfield(handles, 'vof_plots'); end
if isfield(handles, 'best_vof_plot'); delete(handles.best_vof_plot); handles = rmfield(handles, 'best_vof_plot'); end
if isfield(handles, 'vof_pixel_plots'); delete(handles.vof_pixel_plots); handles = rmfield(handles, 'vof_pixel_plots'); end
current_slice = round(get(handles.slider1, 'Value'));
current_timepoint = round(get(handles.slider8, 'Value'));
cla(handles.axes1, 'reset')
handles.vof_conc_image = imshow(zeros(handles.img_size(1), handles.img_size(2)), clim_im, 'colormap', gray, 'Parent', handles.axes1);
set(handles.vof_conc_image, 'CData', (squeeze(handles.dsc_data_c(:,:,current_slice,current_timepoint))) )
cla(handles.axes2, 'reset')
guidata(hObject, handles)
initialise_dynamic_vof_plot(hObject)
initialise_roi_placement(hObject)
end

function saveVOF(hObject, ~, ~)
handles = guidata(hObject);
[file, path] = uiputfile( '*.*','Save VOF', fullfile(handles.file_folder, strcat('vof_',handles.file_name)));
if file == 0
    report('VOF not saved.', handles, 'r')
else
    path_to_vof = fullfile(path,file);
    handles.path_to_vof = path_to_vof;
    guidata(hObject, handles)
    save_data.header_info = [];
    save_data.this_folder = path_to_vof; save_data.this_format = 0;
    save_data.this_name = 0; save_data.data_to_save = handles.best_vof; save_this_file(save_data)
end
end

function mouseClick(hObject, ~, ~)
handles = guidata(hObject);
if ~ishandle(handles.vof_conc_image)||~isfield(handles, 'vof_candidates'); return; end
% disp("Clicked")
C = get(handles.axes1, 'CurrentPoint');
x = round(C(1,1)); %these are flipped
y = round(C(1,2));
pixel_x = x;
pixel_y = y;
if isfield(handles, 'ready_to_edit')
    handles.vof_candidates_copy = handles.vof_candidates;
    handles = rmfield(handles, 'ready_to_edit'); %we only want to get a fresh copy of suggested VOFs once.
end

exist_s = false;
this_indx = false;
for c = 1:handles.numvofs %check if points have already been clicked
    if  handles.opt_vof_coordinates(c, 1) == pixel_x &&  handles.opt_vof_coordinates(c, 2) == pixel_y
        exist_s = true;
        this_indx = c;
    end
end

if exist_s && this_indx ~= false %check if the clicked pixel is one of our optimal pixels
    data_on_pixel = get(handles.vof_pixel_plots(this_indx), 'YData'); %check if its plot exists
    if ~isempty(data_on_pixel) %if it does
        set(handles.vof_pixel_plots(this_indx), 'XData', [], 'YData', []) %delete it
        set( handles.vof_plots(this_indx), 'XData', [], 'YData', []) %and delete the corresponding vof curve
         handles.vof_candidates_copy(this_indx, :) = NaN(1, handles.img_size(4)); %remove the curve from the list of optimal ones
        handles.best_vof = nanmean(handles.vof_candidates_copy, 1); %update the mean
        set(handles.best_vof_plot, 'XData', handles.t, 'YData', handles.best_vof) %and update the mean plot
        drawnow     %now
    else %if the plot nolonger exists (already clicked off)
        set(handles.vof_pixel_plots(this_indx), 'XData', handles.opt_vof_coordinates(this_indx, 1), 'YData', handles.opt_vof_coordinates(this_indx, 2)) %restore it
        set( handles.vof_plots(this_indx), 'XData', handles.t, 'YData', handles.vof_candidates(this_indx, :)) %restore its vof curve
         handles.vof_candidates_copy(this_indx, :) = handles.vof_candidates(this_indx, :); %restore the curve to the list of optimal ones
        handles.best_vof = nanmean(handles.vof_candidates_copy, 1); %update the mean
        set(handles.best_vof_plot, 'XData', handles.t, 'YData', handles.best_vof) %and update the mean plot
        drawnow     %now
    end
else %if clicked pixel is not one of our optimal ones
end % do nothing
handles.vof_area = trapz(handles.t, handles.best_vof);
guidata(hObject, handles)
end

function Exit(hObject, ~, ~)
handles = guidata(hObject);
try
handles = rmfield(handles, {'vof_candidates', 'dsc_data_c', 'vof_plots', 'vof_pixel_plots', 'opt_vof_coordinates', 'roi_slices', 'rois'});
catch
end
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
report('VOF selection closed.', handles)
set([handles.pushbutton24, handles.pushbutton25,handles.pushbutton26, handles.pushbutton27, handles.pushbutton30], 'Visible', 'off')
guidata(hObject, handles)
end