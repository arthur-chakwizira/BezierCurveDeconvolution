function [mean_aif_c, aif_area] = find_aif_auto_gui(handles)
%       This function accepts the handles structure and returns a mean AIF
%       and the area under that mean AIF. It displays a DSC-MRI
%       concentration image on the GUI window and lets the user place a
%       ROI. The function suggests AIFs in the drawn ROI and does not allow
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
aif_handles = handles.aif_handles;
aif_slice = aif_handles.aif_slice;
aif_time_point = aif_handles.aif_time_point;
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
handles.data_is_for_aif = true;
clim_im = [0 100];
handles.aif_conc_image = imshow(zeros(img_size(1), img_size(2)), clim_im, 'colormap', gray, 'Parent', handles.axes1);
handles.aif_conc_data = dsc_data_c;
set(handles.aif_conc_image, 'CData', squeeze(handles.aif_conc_data(:,:,aif_slice,aif_time_point))  )
title(handles.axes1, ['Concentration: Slice ' num2str(aif_slice) ' : Time point ' num2str(aif_time_point)])
set(handles.slider1, 'Min', 1, 'Max', handles.img_size(3), ...
                         'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', aif_slice)
set(handles.slider8, 'Visible', 'on')%horizontal scrolling on axes1
set(handles.slider8, 'Min', 1, 'Max', handles.img_size(4), ...
                         'SliderStep', [1, 1]/(handles.img_size(4) - 1), 'Value', aif_time_point)
clim_max = double(max(handles.aif_conc_data(isfinite(handles.aif_conc_data)))) + 1E-10;
if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/990; end
if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end %$$$$$$$
set(handles.slider6, 'Visible', 'on') %color setting
set(handles.slider6, 'Min', 10E-10, 'Max', clim_max, ...
                         'SliderStep', [step_small, step_big], 'Value', clim_im(2))
set(handles.slider10, 'Visible', 'on') %color setting
set(handles.slider10, 'Min', 0, 'Max', clim_max-1E-10, ...
                         'SliderStep', [step_small, step_big], 'Value', 0)                        
guidata(hObject, handles)

    wrap_text = 'Begin to draw freehand ROI on the image for AIF selection. Use slider to change the displayed slice.';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
for rep = 1:128*128
    if rep==1; blink_this(handles.edit7, 'b'); end
    h_roi = drawfreehand(handles.axes1, 'Color', 'yellow', 'LineWidth', 0.5);
    aif_roi = h_roi.createMask; 
    aif_mask = false(img_size(1:3));
    aif_slice = round(get(handles.slider1, 'Value'));
    
    aif_mask(:,:,aif_slice) = aif_roi.*mask(:,:,aif_slice);
 
tmp_dsc_data_c = reshape(dsc_data_c,[prod(img_size(1:3)) img_size(4)]); 
tmp_mask = repmat(aif_mask(:),[1 img_size(4)]);
aif_data = reshape(tmp_dsc_data_c(tmp_mask),[sum(aif_mask(:)) img_size(4)]);
aif_data_s = reshape(dsc_data(tmp_mask),[sum(aif_mask(:)) img_size(4)]);

if sum(aif_mask(:)) < 1
       wrap_text = 'Current AIF selection method does not support point-ROIs. Please draw larger ROI (that includes non-background pixels).';
       set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
       for rr = 1:5
           h_roi.Color = 'black';
           pause(0.15)
           h_roi.Color = 'red';
           pause(0.15)
       end
       delete(h_roi)
else

for ii = 1:sum(aif_mask(:))    
    [aif_max, tmax] = max(aif_data(ii,:)); %maximum of AIF concentration and time point at which it occurs
    aif_baseline = mean(aif_data_s(ii, baseline_idx)); %baseline of AIF signal
    aif_tail = mean(aif_data_s(ii, tail_idx)); %tail of the AIF signal
    t20percent = find(aif_data(ii,:) > 0.2*aif_max); %find where AIF concentration is above 20% of its maximum (that is, locate peak)
    if isempty(t20percent) %happens if curve is just a bunch of zeros
        aif_data(ii, :) = NaN; %discriminate against this
%         disp('no peak')
    else
        num_peaks_check = diff(t20percent); %checking for curves with multiple distinct maxima
        if ~all(num_peaks_check == 1) %if there is a jump of more than 1 in the peak indices, the curve has more than one peak
            aif_data(ii, :) = NaN; %discriminate
%             disp('2 peaks')
        else %if there is just one peak (it may still have small subpeaks that are not detected by the check above)
%             disp('1 peak with properties: \n')
            t20_low = t20percent(1); %extract beginning
            t20_up = t20percent(end); %and end of peak
            rising_end = diff(aif_data(ii, t20_low:tmax)); %check slope on rising side
            falling_end = diff(aif_data(ii, tmax:t20_up));  %check slope on falling side
            if ~all(rising_end > 0); aif_data(ii, :) = NaN; end %if not monotonically increasing, discriminate
            if ~all(falling_end < 0); aif_data(ii, :) = NaN; end %if not monotonically decreasing, discriminate
        end
    end
%     if aif_baseline < 0.6*max_baseline; aif_c(ii, :) = 0; end %if AIF signal baseline is lower than 60% of maximum in ROI, discriminate
%     if aif_tail < 0.8*aif_baseline; aif_c(ii, :) = 0; end %if AIF signal tail is lower than 80% of original baseline, discriminate
    if min(aif_data_s(ii,:)) > 0.2*aif_baseline; aif_data(ii, :) = NaN; end %if AIF signal drop is lower than 80% of baseline, discriminate
    if any(aif_data(ii, :) < -10); aif_data(ii, :) = NaN; end
end

aif_data_c = aif_data;
ymax    = zeros(sum(aif_mask(:)),1);    
fm      = zeros(sum(aif_mask(:)),1);  
noise   = zeros(sum(aif_mask(:)),1);
pointy  = zeros(sum(aif_mask(:)),1);
for i = 1:sum(aif_mask(:))
    if all(isnan(aif_data_c(i, :))) %if AIF has previously been discriminated against
        ymax(i) = 0; fm(i) = t(end); noise(i) = 1E6; pointy(i) = 1E6; %discriminate
    else    %otherwise    
    [ymax(i), fm(i)] = max(aif_data_c(i,:)); 
%     fm(i) = sum(t.*aif_data_c(i,:))/sum(aif_data_c(i,:));
    %fm(i) = sum(t(fit_idx).*aif_data_c(i,fit_idx))/sum(aif_data_c(i,fit_idx));
    % noisy
    noise(i) = std(aif_data(i,baseline_idx)) + std(aif_data(i,tail_idx));
    % pointyness
    pointy(i) = sum(abs(diff(aif_data(i,:))));    
    end
end

numaifs = 10; %number of curves to work with
if numaifs > sum(aif_mask(:)); numaifs = sum(aif_mask(:)); end %if greater than number of selected pixels, set to number of selected pixels
[~,sortIndex] = sort(ymax(:),'descend'); %sort in descending order
ymax_index = sortIndex(1:numaifs); %get highest 10

[~,sortIndex] = sort(noise(:),'ascend');  %sort in ascending noise
noise_index = sortIndex(1:numaifs); % get lowest 10

[~,sortIndex] = sort(pointy(:),'ascend'); %sort in scending pointyness
pointy_index = sortIndex(1:numaifs); % get lowest 10

[~,sortIndex] = sort(fm(:),'ascend'); % sort in scending fm
fm_index = sortIndex(1:numaifs); % get lowest 10


% find intersects
ymaxfm_index = intersect(intersect(intersect(ymax_index, fm_index), noise_index), pointy_index);%initially demand satisfaction of all requirements
if (isempty(ymaxfm_index)||length(ymaxfm_index)<5); ymaxfm_index =intersect(intersect(ymax_index, fm_index), noise_index); end %if failed, remove low pointyness requirement
if (isempty(ymaxfm_index)||length(ymaxfm_index)<5); ymaxfm_index = intersect(ymax_index, fm_index); end %if still no candidates, remove low fm requirement

numaifs = length(ymaxfm_index); %check how many AIFs we have left 
if numaifs < 5 %if none
ymaxfm_index = fm_index;  %use only those with max value, that is, remove low noise requirement
numaifs = length(ymaxfm_index); %check how many we have left now
end
if numaifs > 5; numaifs = 5; ymaxfm_index = ymaxfm_index(1:5); end %if more than 5, pick only the first 5

h_roi.Color = 'green';
mean_aif_c = nanmean(aif_data_c(ymaxfm_index,:),1);
plot(t, mean_aif_c, 'k-', 'LineWidth', 2, 'Parent', handles.axes2);
hold(handles.axes2, 'on')
drawnow
for ai_f = 1:numaifs
   plot(t, aif_data_c(ymaxfm_index(ai_f), :),'HandleVisibility', 'off','LineWidth', 0.8, 'Parent', handles.axes2) 
end
%     xlim(handles.axes2, [0 img_size(4)])
%     ylim(handles.axes2, [-10 150])
%     set(handles.axes2, 'YTick', -10:10:150)
    xlabel(handles.axes2, 'time [s]')
    ylabel(handles.axes2, 'Concentration')
    title(handles.axes2, 'AIFs in selected ROI')
    legend(handles.axes2, 'Mean')
  aif_area = trapz(t,mean_aif_c);  
 guidata(hObject, handles)

 %-----------------------------------------------------------------------
   wrap_text = ['Press S/Enter to save mean AIF. Press C to save and proceed with current AIF.'...
       'Press Esc to terminate selection. Click on image to redo selection.'];
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    blink_this(handles.edit7, 'b')
 user_pressed = waitforbuttonpress;
 this_key = double(get(handles.figure1, 'CurrentCharacter'));
 if user_pressed
            if this_key == 115 || this_key == 13 %s or enter
                [file,path] = uiputfile(file_filter);
                if file == 0
                    wrap_text = 'AIF not saved. Selection still active.';
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                    pause(2)
                    wrap_text = 'Begin to draw freehand ROI on the image for AIF selection. Use slider to change the displayed slice.';
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
                    if ishandle(h_roi); delete(h_roi); end
                    cla(handles.axes2, 'reset')
                    xlim(handles.axes2, [0 img_size(4)])
%                     ylim(handles.axes2, [-10 150])
                    xlabel(handles.axes2, 'time [s]')
                    ylabel(handles.axes2, 'Concentration')
                    title(handles.axes2, 'AIFs in selected ROI')
                else
                path_to_aif = fullfile(path,file);
                try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
                save_data.this_folder = path_to_aif; save_data.this_format = 0;
                save_data.this_name = 0; save_data.data_to_save = mean_aif_c; save_this_file(save_data)
                end
            end
            if this_key == 27
                mean_aif_c = false; aif_area = false;
                handles.data_is_for_aif = false;
                set(handles.pushbutton14, 'Visible', 'off')
                guidata(hObject, handles)
                return
            end %escape
            if this_key == 99 %c
                aif_area = trapz(t, mean_aif_c);
                [file,path] = uiputfile(file_filter);
                if file == 0
                    wrap_text = 'AIF not saved. Analysis is proceeding...';
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                    pause(3)
                else
                path_to_aif = fullfile(path,file);
                try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
                save_data.this_folder = path_to_aif; save_data.this_format = 0;
                save_data.this_name = 0; save_data.data_to_save = mean_aif_c; save_this_file(save_data)
                end
                handles.data_is_for_aif = false;
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
    title(handles.axes2, 'AIFs in selected ROI')
 end%c

    end
end
handles.data_is_for_aif = false;
guidata(hObject, handles)
    
end