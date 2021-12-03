
function [mean_vof_c, vof_area] = find_vof_auto_gui(handles)
%       This function accepts the handles structure and returns a mean VOF
%       and the area under that mean VOF. It displays a DSC-MRI
%       concentration image on the GUI window and lets the user place a
%       ROI. The function suggests VOFs in the drawn ROI and does not allow
%       the user to edit the suggestion. ROI placement can be repeated a
%       large number of times (16384).
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dsc_data = handles.dsc_data;
mask = handles.mask;
te = handles.te; %s
tr =  handles.tr; %s 
img_size = handles.img_size;
baseline_idx = handles.baseline_index;
tail_idx = handles.tail_index;
vof_handles = handles.vof_handles;
vof_slice = vof_handles.vof_slice;
vof_time_point = vof_handles.vof_time_point;
t = 0:tr:(img_size(4)-1)*tr;
hObject = handles.hObject;
file_filter = handles.file_filter;

dsc_data_c = dsc_data.*0;
for ro_w = 1:img_size(1)
    for co_l = 1:img_size(2)
        for sl_c = 1:img_size(3)
            if mask(ro_w, co_l, sl_c)
                pixel_signal = squeeze(dsc_data(ro_w, co_l, sl_c, :));
                pixel_s0 = mean(pixel_signal(baseline_idx));
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
if ishandle(handles.wait_for_data); delete(handles.wait_for_data); end
handles.data_is_for_vof = true;
clim_im = [0 100];
handles.vof_conc_image = imshow(zeros(img_size(1), img_size(2)), clim_im, 'colormap', gray, 'Parent', handles.axes1);
handles.vof_conc_data = dsc_data_c;
set(handles.vof_conc_image, 'CData', squeeze(handles.vof_conc_data(:,:,vof_slice,vof_time_point))  )
title(handles.axes1, ['Concentration: Slice ' num2str(vof_slice) ' : Time point ' num2str(vof_time_point)])
set(handles.slider1, 'Min', 1, 'Max', handles.img_size(3), ...
                         'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', vof_slice)
set(handles.slider8, 'Visible', 'on')%horizontal scrolling on axes1
set(handles.slider8, 'Min', 1, 'Max', handles.img_size(4), ...
                         'SliderStep', [1, 1]/(handles.img_size(4) - 1), 'Value', vof_time_point)
clim_max = double(max(handles.vof_conc_data(isfinite(handles.vof_conc_data)))) + 1E-10;
if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/990; end
if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end %$$$$$$$
set(handles.slider6, 'Visible', 'on') %color setting
set(handles.slider6, 'Min', 10E-10, 'Max', clim_max, ...
                         'SliderStep', [step_small, step_big], 'Value', clim_im(2))
set(handles.slider10, 'Visible', 'on') %color setting
set(handles.slider10, 'Min', 0, 'Max', clim_max-1E-10, ...
                         'SliderStep', [step_small, step_big], 'Value', 0)                        
guidata(hObject, handles)

    wrap_text = 'Begin to draw freehand ROI on the image for VOF selection. Use slider to change the displayed slice.';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
for rep = 1:128*128
    if rep==1; blink_this(handles.edit7, 'b'); end
    h_roi = drawfreehand(handles.axes1, 'Color', 'yellow', 'LineWidth', 0.5);
    vof_roi = h_roi.createMask; 
    vof_mask = false(img_size(1:3));
    vof_slice = round(get(handles.slider1, 'Value'));
    
    vof_mask(:,:,vof_slice) = vof_roi.*mask(:,:,vof_slice);
 
tmp_dsc_data_c = reshape(dsc_data_c,[prod(img_size(1:3)) img_size(4)]); 
tmp_mask = repmat(vof_mask(:),[1 img_size(4)]);
vof_data = reshape(tmp_dsc_data_c(tmp_mask),[sum(vof_mask(:)) img_size(4)]);
vof_data_s = reshape(dsc_data(tmp_mask),[sum(vof_mask(:)) img_size(4)]);

if sum(vof_mask(:)) < 1
       wrap_text = 'Current VOF selection method does not support point-ROIs. Please draw larger ROI (that includes non-background pixels).';
       set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
       for rr = 1:5
           h_roi.Color = 'black';
           pause(0.15)
           h_roi.Color = 'red';
           pause(0.15)
       end
       delete(h_roi)
else

for ii = 1:sum(vof_mask(:))    
    [vof_max, tmax] = max(vof_data(ii,:)); %maximum of VOF concentration and time point at which it occurs
    vof_baseline = mean(vof_data_s(ii, baseline_idx)); %baseline of VOF signal
    vof_tail = mean(vof_data_s(ii, tail_idx)); %tail of the VOF signal
    t20percent = find(vof_data(ii,:) > 0.2*vof_max); %find where VOF concentration is above 20% of its maximum (that is, locate peak)
    if isempty(t20percent) %happens if curve is just a bunch of zeros
        vof_data(ii, :) = NaN; %discriminate against this
%         disp('no peak')
    else
        num_peaks_check = diff(t20percent); %checking for curves with multiple distinct maxima
        if ~all(num_peaks_check == 1) %if there is a jump of more than 1 in the peak indices, the curve has more than one peak
            vof_data(ii, :) = NaN; %discriminate
%             disp('2 peaks')
        else %if there is just one peak (it may still have small subpeaks that are not detected by the check above)
%             disp('1 peak with properties: \n')
            t20_low = t20percent(1); %extract beginning
            t20_up = t20percent(end); %and end of peak
            rising_end = diff(vof_data(ii, t20_low:tmax)); %check slope on rising side
            falling_end = diff(vof_data(ii, tmax:t20_up));  %check slope on falling side
            if ~all(rising_end > 0); vof_data(ii, :) = NaN; end %if not monotonically increasing, discriminate
            if ~all(falling_end < 0); vof_data(ii, :) = NaN; end %if not monotonically decreasing, discriminate
        end
    end
%     if vof_baseline < 0.6*max_baseline; vof_c(ii, :) = 0; end %if VOF signal baseline is lower than 60% of maximum in ROI, discriminate
%     if vof_tail < 0.8*vof_baseline; vof_c(ii, :) = 0; end %if VOF signal tail is lower than 80% of original baseline, discriminate
    if min(vof_data_s(ii,:)) > 0.3*vof_baseline; vof_data(ii, :) = NaN; end %if VOF signal drop is lower than 80% of baseline, discriminate
    if any(vof_data(ii, :) < -10); vof_data(ii, :) = NaN; end
end

vof_data_c = vof_data;
ymax    = zeros(sum(vof_mask(:)),1);    
fm      = zeros(sum(vof_mask(:)),1);  
noise   = zeros(sum(vof_mask(:)),1);
pointy  = zeros(sum(vof_mask(:)),1);
for i = 1:sum(vof_mask(:))
    if all(isnan(vof_data_c(i, :))) %if VOF has previously been discriminated against
        ymax(i) = 0; fm(i) = t(1); noise(i) = 1E6; pointy(i) = 1E6; %discriminate
    else    %otherwise    
    [ymax(i), fm(i)] = max(vof_data_c(i,:)); 
%     fm(i) = sum(t.*vof_data_c(i,:))/sum(vof_data_c(i,:));
    %fm(i) = sum(t(fit_idx).*vof_data_c(i,fit_idx))/sum(vof_data_c(i,fit_idx));
    % noisy
    noise(i) = std(vof_data(i,baseline_idx)) + std(vof_data(i,tail_idx));
    % pointyness
    pointy(i) = sum(abs(diff(vof_data(i,:))));    
    end
end

numvofs = 10; %number of curves to work with
if numvofs > sum(vof_mask(:)); numvofs = sum(vof_mask(:)); end %if greater than number of selected pixels, set to number of selected pixels
[~,sortIndex] = sort(ymax(:),'descend'); %sort in descending order
ymax_index = sortIndex(1:numvofs); %get highest 10

[~,sortIndex] = sort(noise(:),'ascend');  %sort in ascending noise
noise_index = sortIndex(1:numvofs); % get lowest 10

[~,sortIndex] = sort(pointy(:),'ascend'); %sort in scending pointyness
pointy_index = sortIndex(1:numvofs); % get lowest 10

[~,sortIndex] = sort(fm(:),'descend'); % sort in scending fm
fm_index = sortIndex(1:numvofs); % get lowest 10


% find intersects
ymaxfm_index = intersect(intersect(intersect(ymax_index, fm_index), noise_index), pointy_index);%initially demand satisfaction of all requirements
if (isempty(ymaxfm_index)||length(ymaxfm_index)<5); ymaxfm_index =intersect(intersect(ymax_index, fm_index), noise_index); end %if failed, remove low pointyness requirement
if (isempty(ymaxfm_index)||length(ymaxfm_index)<5); ymaxfm_index = intersect(ymax_index, fm_index); end %if still no candidates, remove low fm requirement

numvofs = length(ymaxfm_index); %check how many VOFs we have left 
if numvofs < 5 %if none
ymaxfm_index = fm_index;  %use only those with max value, that is, remove low noise requirement
numvofs = length(ymaxfm_index); %check how many we have left now
end
if numvofs > 5; numvofs = 5; ymaxfm_index = ymaxfm_index(1:5); end %if more than 5, pick only the first 5

h_roi.Color = 'green';
mean_vof_c = nanmean(vof_data_c(ymaxfm_index,:),1);
plot(t, mean_vof_c, 'k-', 'LineWidth', 2, 'Parent', handles.axes2);
hold(handles.axes2, 'on')
for ai_f = 1:numvofs
   plot(t, vof_data_c(ymaxfm_index(ai_f), :),'HandleVisibility', 'off','LineWidth', 0.8, 'Parent', handles.axes2) 
end
%     xlim(handles.axes2, [0 img_size(4)])
%     ylim(handles.axes2, [-10 150])
%     set(handles.axes2, 'YTick', -10:10:150)
    xlabel(handles.axes2, 'time [s]')
    ylabel(handles.axes2, 'Concentration')
    title(handles.axes2, 'VOFs in selected ROI')
    legend(handles.axes2, 'Mean')
  vof_area = trapz(t,mean_vof_c);  
 guidata(hObject, handles)

 %-----------------------------------------------------------------------
   wrap_text = ['Press S/Enter to save mean VOF. Press C to save and proceed with current VOF.'...
       'Press Esc to terminate selection. Click on image to redo selection.'];
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    blink_this(handles.edit7, 'b')
 user_pressed = waitforbuttonpress;
 this_key = double(get(handles.figure1, 'CurrentCharacter'));
 if user_pressed
            if this_key == 115 || this_key == 13 %s or enter
                [file,path] = uiputfile(file_filter);
                if file == 0
                    wrap_text = 'VOF not saved. Selection still active.';
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                    pause(2)
                    wrap_text = 'Begin to draw freehand ROI on the image for VOF selection. Use slider to change the displayed slice.';
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
                    if ishandle(h_roi); delete(h_roi); end
                    cla(handles.axes2, 'reset')
                    xlim(handles.axes2, [0 img_size(4)])
%                     ylim(handles.axes2, [-10 150])
                    xlabel(handles.axes2, 'time [s]')
                    ylabel(handles.axes2, 'Concentration')
                    title(handles.axes2, 'VOFs in selected ROI')
                else
                path_to_vof = fullfile(path,file);
                try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
                save_data.this_folder = path_to_vof; save_data.this_format = 0;
                save_data.this_name = 0; save_data.data_to_save = mean_vof_c; save_this_file(save_data)
                end
            end
            if this_key == 27
                mean_vof_c = false; vof_area = false;
                handles.data_is_for_vof = false;
                set(handles.pushbutton14, 'Visible', 'off')
                guidata(hObject, handles)
                return
            end %escape
            if this_key == 99 %c
                vof_area = trapz(t, mean_vof_c);
                [file,path] = uiputfile(file_filter);
                if file == 0
                    wrap_text = 'VOF not saved. Analysis is proceeding...';
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                    pause(3)
                else
                path_to_vof = fullfile(path,file);
                try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
                save_data.this_folder = path_to_vof; save_data.this_format = 0;
                save_data.this_name = 0; save_data.data_to_save = mean_vof_c; save_this_file(save_data)
                end
                handles.data_is_for_vof = false;
                set(handles.pushbutton14, 'Visible', 'off')
                guidata(hObject, handles)
                return
            end
 else
    if ishandle(h_roi); delete(h_roi); end
    cla(handles.axes2, 'reset')
        xlim(handles.axes2, [0 img_size(4)])
%     ylim(handles.axes2, [-10 150])
    xlabel(handles.axes2, 'time [s]')
    ylabel(handles.axes2, 'Concentration')
    title(handles.axes2, 'VOFs in selected ROI')
 end%c

    end
end
handles.data_is_for_vof = false;
guidata(hObject, handles)
    
end