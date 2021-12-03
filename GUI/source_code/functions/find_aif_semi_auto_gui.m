function [mean_aif_c, aif_area] = find_aif_semi_auto_gui(handles)
%       This function accepts the handles structure and returns a mean AIF
%       and the area under that mean AIF. It displays a DSC-MRI
%       concentration image on the GUI window and lets the user place a
%       ROI. The function suggests AIFs in the drawn ROI and allows the
%       user to edit the suggestion.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hObject = handles.hObject;
file_filter = handles.file_filter;
dsc_data = handles.dsc_data;
mask = handles.mask;
te = handles.te; %s
tr = handles.tr; %s 
img_size = handles.img_size;
baseline_idx = handles.baseline_index;
tail_idx = handles.tail_index;
aif_handles = handles.aif_handles;
aif_slice = aif_handles.aif_slice;
aif_time_point = aif_handles.aif_time_point;
t = 0:tr:(img_size(4)-1)*tr;
%convert dsc_data to concentration
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
if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end 
set(handles.slider6, 'Visible', 'on') %color setting
set(handles.slider6, 'Min', 10E-10, 'Max', clim_max, ...
                         'SliderStep', [step_small, step_big], 'Value', clim_im(2))
set(handles.slider10, 'Visible', 'on') %color setting
set(handles.slider10, 'Min', 0, 'Max', clim_max-1E-10, ...
                         'SliderStep', [step_small, step_big], 'Value', 0)                     
guidata(hObject, handles)

wrap_text = ['Click on image to begin drawing of freehand ROI for AIF selection.  '...  
    'Program will automatically suggest AIFs in selected ROI.   Use slider to change displayed slice.' ];
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
blink_this(handles.edit7, 'b')
stop_now = waitforbuttonpress;
if stop_now ~= 0
%     mean_aif_c = false; aif_area = false; return;
    dont_key = double(get(handles.figure1,'CurrentCharacter'));
    if dont_key == 27; mean_aif_c = false; aif_area = false; return; end
end

for rep = 1:128^4
%     try
    saving_now = false;
    if rep >1
        stop_now = waitforbuttonpress;
        if ~stop_now %if mouse clicked, continue
        else
            just_stop = double(get(handles.figure1,'CurrentCharacter'));
            if just_stop == 27; mean_aif_c = false; aif_area = false; return; end
            if just_stop == 99 %c
                mean_aif_c = opt_mean_aif_c; aif_area = trapz(t, mean_aif_c);
                if ishandle(h_roi); delete(h_roi); end
                    for plt = 1:numaifs
                        if ishandle(opt_roi_plots(plt)); delete(opt_roi_plots(plt)); end
                        if ishandle(opt_aif_plots(plt)); delete(opt_aif_plots(plt)); end    
                    end
                    hold(handles.axes1, 'off')
                    hold(handles.axes2, 'off')
                    [file,path] = uiputfile(file_filter);
                    if file == 0
                        wrap_text = 'AIF not saved. Analysis proceeding...';
                        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                        pause(3)
                   else   
                        path_to_aif = fullfile(path,file);
                        try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
                        save_data.this_folder = path_to_aif; save_data.this_format = 0;
                        save_data.this_name = 0; save_data.data_to_save = mean_aif_c; save_this_file(save_data)
                   end
                   return
            end
        end
     end
    h_roi = drawfreehand(handles.axes1, 'Color', 'yellow' , 'LineWidth', 1);
    h_roi.FaceAlpha = 0.1; %ROI transparency
    h_roi.FaceSelectable = false; %you cannot drag this ROI
    
    %extract coordinates of all pixels in the drawn ROI
    pixel_xcoords = NaN(img_size(1)^2,1);
    pixel_ycoords = NaN(img_size(2)^2,1);
    this_pos = 0;
    for x = 1:img_size(1)
        for y = 1:img_size(2)
            if inROI(h_roi, x,y)
                this_pos = this_pos + 1;
                pixel_xcoords(this_pos) = x;
                pixel_ycoords(this_pos) = y;
            end
        end
    end
    %delete nan values
    pixel_xcoords(isnan(pixel_xcoords)) = [];
    pixel_ycoords(isnan(pixel_ycoords)) = [];
    %forbid point ROIs
    if (isempty(pixel_xcoords) || length(pixel_xcoords) < 5) %@@@@@
       wrap_text = 'Current AIF selection method does not support point-ROIs. Please draw larger ROI.';
       set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
       for rr = 1:5
           h_roi.Color = 'black';
           pause(0.15)
           h_roi.Color = 'red';
           pause(0.15)
       end
       delete(h_roi)
    else
    aif_c = NaN(length(pixel_xcoords), img_size(4));
    aif_s = NaN(length(pixel_xcoords), img_size(4));%@@
    aif_slice = round(get(handles.slider1, 'Value'));
    for k = 1:length(pixel_xcoords)
        aif_s(k,:) = squeeze(dsc_data(pixel_ycoords(k), pixel_xcoords(k), aif_slice, :));
        aif_c(k,:) = squeeze(dsc_data_c(pixel_ycoords(k), pixel_xcoords(k), aif_slice, :));        
    end

ymax    = zeros(length(pixel_xcoords),1);    
fm      = zeros(length(pixel_xcoords),1);  
noise   = zeros(length(pixel_xcoords),1);
pointy  = zeros(length(pixel_xcoords),1);

max_baseline = max(mean(aif_s(:, baseline_idx), 2));

for ii = 1:length(pixel_xcoords) %loop through all pixels in ROI
    [aif_max, tmax] = max(aif_c(ii,:)); %maximum of AIF concentration and time point at which it occurs
    aif_baseline = mean(aif_s(ii, baseline_idx)); %baseline of AIF signal
    aif_tail = mean(aif_s(ii, tail_idx)); %tail of the AIF signal
    t20percent = find(aif_c(ii,:) > 0.2*aif_max); %find where AIF concentration is above 20% of its maximum (that is, locate peak)
    if isempty(t20percent) %happens if curve is just a bunch of zeros
        aif_c(ii, :) = NaN; %discriminate against this
%         disp('no peak')
    else
        num_peaks_check = diff(t20percent); %checking for curves with multiple distinct maxima
        if ~all(num_peaks_check == 1) %if there is a jump of more than 1 in the peak indices, the curve has more than one peak
            aif_c(ii, :) = NaN; %discriminate
%             disp('2 peaks')
        else %if there is just one peak (it may still have small subpeaks that are not detected by the check above)
%             disp('1 peak with properties: \n')
            t20_low = t20percent(1); %extract beginning
            t20_up = t20percent(end); %and end of peak
            rising_end = diff(aif_c(ii, t20_low:tmax)); %check slope on rising side
            falling_end = diff(aif_c(ii, tmax:t20_up));  %check slope on falling side
            if ~all(rising_end > 0); aif_c(ii, :) = NaN; end %if not monotonically increasing, discriminate
            if ~all(falling_end < 0); aif_c(ii, :) = NaN; end %if not monotonically decreasing, discriminate
        end
    end
%     if aif_baseline < 0.6*max_baseline; aif_c(ii, :) = 0; end %if AIF signal baseline is lower than 60% of maximum in ROI, discriminate
%     if aif_tail < 0.8*aif_baseline; aif_c(ii, :) = 0; end %if AIF signal tail is lower than 80% of original baseline, discriminate
    if min(aif_s(ii,:)) > 0.2*aif_baseline; aif_c(ii, :) = NaN; end %if AIF signal drop is lower than 80% of baseline, discriminate
    if any(aif_c(ii, :) < -10); aif_c(ii,:) = NaN; end
end

aif_data_c = aif_c;
for i = 1:length(pixel_xcoords)
    if all(isnan(aif_data_c(i, :))) %if AIF has previously been discriminated against
        ymax(i) = 0; fm(i) = t(end); noise(i) = 1E6; pointy(i) = 1E6; %discriminate
    else    %otherwise
    [ymax(i), fm(i)] = max(aif_data_c(i,:)); %extract amplitude and time of occurence
%     fm(i) = sum(t.*aif_data_c(i,:))/sum(aif_data_c(i,:));
    %fm(i) = sum(t(fit_idx).*aif_data_c(i,fit_idx))/sum(aif_data_c(i,fit_idx));
    noise(i) = std(aif_data_c(i,baseline_idx)) + std(aif_data_c(i,tail_idx)); %compute noise level (baseline + tail)
    pointy(i) = sum(abs(diff(aif_data_c(i,:)))); %compute pointyness
    end
end

numaifs = 10; %number of curves to work with
if numaifs > length(pixel_xcoords); numaifs = length(pixel_xcoords); end %if greater than number of selected pixels, set to number of selected pixels
[~,sortIndex] = sort(ymax(:),'descend'); %sort in descending order
ymax_index = sortIndex(1:numaifs); % get highest 10

[~,sortIndex] = sort(noise(:),'ascend'); %sort in ascending noise
noise_index = sortIndex(1:numaifs); % get lowest 10

[~,sortIndex] = sort(pointy(:),'ascend'); % sort in scending pointyness
pointy_index = sortIndex(1:numaifs); %get lowest 10

[~,sortIndex] = sort(fm(:),'ascend'); %sort in ascending fm
fm_index = sortIndex(1:numaifs); % get lowest 10


% find intersects
ymaxfm_index = intersect(intersect(intersect(ymax_index, fm_index), noise_index), pointy_index); %initially demand satisfaction of all requirements
if (isempty(ymaxfm_index)||length(ymaxfm_index)<5); ymaxfm_index =intersect(intersect(ymax_index, fm_index), noise_index); end %if failed, remove low pointyness requirement
if (isempty(ymaxfm_index)||length(ymaxfm_index)<5); ymaxfm_index = intersect(ymax_index, fm_index); end %if still no candidates, remove noise requirement

numaifs = length(ymaxfm_index); %check how many AIFs we have left 
if numaifs < 5 %if none
ymaxfm_index = fm_index; %use only those with earliest arrival, that is, remove all other requirements
numaifs = length(ymaxfm_index); %check how many we have left now
end
if numaifs > 5; numaifs = 5; ymaxfm_index = ymaxfm_index(1:5); end %if more than 5, pick only the first 5

    %plot curves deemed optimal
    h_roi.Color = 'green';
    hold(handles.axes1, 'on')
    opt_mean_aif_plot = plot(0, 0, 'k-','LineWidth', 2.2, 'Parent', handles.axes2);
%     xlim(handles.axes2, [0 img_size(4)])
%     ylim(handles.axes2, [-10 150])
%     set(handles.axes2, 'YTick', -10:10:150)
    xlabel(handles.axes2, 'time [s]')
    ylabel(handles.axes2, 'Concentration')
    title(handles.axes2, 'AIF selection')
    legend(handles.axes2, 'Mean')
    hold(handles.axes2, 'on')
    opt_roi_plots = zeros(numaifs, 1); 
    opt_aif_plots = zeros(numaifs,1);
    opt_aif_c = NaN(numaifs, img_size(4));
    opt_coordinates = zeros(numaifs, 2); %x y
%     opt_colors = zeros(numaifs,1);
    color_list = {rgb('RoyalBlue'), rgb('Red'), rgb('Green'), rgb('Orange'), rgb('Cyan')};
    for k = 1:numaifs
        opt_aif_c(k,:) = aif_data_c(ymaxfm_index(k), :);
        opt_aif_plots(k) = plot(t, opt_aif_c(k, :),'color',color_list{k},'LineWidth', 0.8,'HandleVisibility', 'off', 'Parent', handles.axes2);        
        this_color = get(opt_aif_plots(k), 'color');
%         opt_colors(k) = this_color; %we will save these for later
        if all(isnan(opt_aif_c(k, :))); marker_size = 0.01; else; marker_size = 5; end %@@@
        opt_roi_plots(k) = plot(pixel_xcoords(ymaxfm_index(k)), pixel_ycoords(ymaxfm_index(k)), 's',...
            'MarkerFaceColor', this_color, 'MarkerEdgeColor', this_color, 'MarkerSize', marker_size);
        opt_mean_aif_c = nanmean(opt_aif_c(1:k, :), 1);
        set(opt_mean_aif_plot, 'XData', t, 'YData', opt_mean_aif_c)
        drawnow
        opt_coordinates(k,:) = [pixel_xcoords(ymaxfm_index(k)),   pixel_ycoords(ymaxfm_index(k))];
    end

wrap_text = ['Press:  R or Spacebar to redo ROI placement,   S or Enter to save current AIF, '...
    ' C to continue with currently displayed AIF, '...
    ' Esc to terminate AIF selection completely.   Click anywhere on the window to edit AIF selection. '];
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
blink_this(handles.edit7, 'b') 

redraw_roi = waitforbuttonpress; %allow redrawing of roi if user presses any key; continue if mouse click
clear('key_pushed')
key_pushed = double(get(handles.figure1,'CurrentCharacter'));

    if ~redraw_roi

%call ginput to allow deselection of individual pixels.
    n_clicks = 128*128; %maximum number of clicks is this. We don't expect more :-)
opt_aif_c_copy = opt_aif_c;

    for i = 1:n_clicks
        wrap_text = ['Click once to activate cursor.  Click on pixel to deselect/reselect it.    '...  
            'Press: R or Spacebar to redraw ROI,   S or Enter to save AIF,   Esc to exit AIF selection,   C to continue with current AIF.'];
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
        if i == 1; blink_this(handles.edit7, 'b'); end 
        if ishandle(h_roi); delete(h_roi); end
        what_now = waitforbuttonpress; %0 for mouse-click, 1 for key press
        if ~what_now
        else
        clear('the_key')
        the_key = double(get(handles.figure1,'CurrentCharacter'));
            if the_key == 114 || the_key == 32 %r or spacebar
                if ishandle(h_roi); delete(h_roi); end
                for plt = 1:numaifs
                    if ishandle(opt_roi_plots(plt)); delete(opt_roi_plots(plt)); end
                    if ishandle(opt_aif_plots(plt)); delete(opt_aif_plots(plt)); end
                end
                hold(handles.axes1, 'off')
                hold(handles.axes2, 'off')                
%                 wrap_text = ['Click on image to begin drawing of freehand ROI for AIF selection.  '...  
%                 'Program will automatically suggest AIFs in selected ROI.   Use slider to change displayed slice.' ];
%                 set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
%                 blink_this(handles.edit7, 'b')           
                break
            end
            if the_key == 115 || the_key == 13 %s or enter
                [file,path] = uiputfile(file_filter);
                if file == 0
                    wrap_text = 'AIF not saved. Selection is still active.';
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                    pause(2) 
                else
                path_to_aif = fullfile(path,file);
                try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
                save_data.this_folder = path_to_aif; save_data.this_format = 0;
                save_data.this_name = 0; save_data.data_to_save = mean_aif_c; save_this_file(save_data)
                saving_now = true;
                break
                end
            end
            if the_key == 27 %escape
                    mean_aif_c = false;
                    aif_area = false;
                    return
            end
            if the_key == 99 %c
                if ishandle(h_roi); delete(h_roi); end
                for plt = 1:numaifs
                    if ishandle(opt_roi_plots(plt)); delete(opt_roi_plots(plt)); end
                    if ishandle(opt_aif_plots(plt)); delete(opt_aif_plots(plt)); end    
                end
                hold(handles.axes1, 'off')
                hold(handles.axes2, 'off')
                [file,path] = uiputfile(file_filter);
                if file == 0
                    wrap_text = 'AIF not saved. Analysis proceeding...';
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                    pause(3)
                else   
                    path_to_aif = fullfile(path,file);
                    try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
                    save_data.this_folder = path_to_aif; save_data.this_format = 0;
                    save_data.this_name = 0; save_data.data_to_save = mean_aif_c; save_this_file(save_data)
                end
                return
            end   
        end

        [x, y] = the_ginput(1);
        pixel_x = round(x);
        pixel_y = round(y);
        exist_s = false;
        this_indx = false;
        for c = 1:numaifs %check if points have already been clicked
            if opt_coordinates(c, 1) == pixel_x & opt_coordinates(c, 2) == pixel_y
                exist_s = true;
                this_indx = c;
            end
        end
        
        if exist_s & this_indx ~= false %check if the clicked pixel is one of our optimal pixels
            data_on_pixel = get(opt_roi_plots(this_indx), 'YData'); %check if its plot exists
            if ~isempty(data_on_pixel) %if it does
                set(opt_roi_plots(this_indx), 'XData', [], 'YData', []) %delete it
                set(opt_aif_plots(this_indx), 'XData', [], 'YData', []) %and delete the corresponding aif curve
                opt_aif_c_copy(this_indx, :) = NaN(1, img_size(4)); %remove the curve from the list of optimal ones
                opt_mean_aif_c = nanmean(opt_aif_c_copy, 1); %update the mean
                set(opt_mean_aif_plot, 'XData', t, 'YData', opt_mean_aif_c) %and update the mean plot
                drawnow     %now
            else %if the plot nolonger exists (already clicked off)
                set(opt_roi_plots(this_indx), 'XData', opt_coordinates(this_indx, 1), 'YData', opt_coordinates(this_indx, 2)) %restore it
                set(opt_aif_plots(this_indx), 'XData', t, 'YData', opt_aif_c(this_indx, :)) %restore its aif curve
                opt_aif_c_copy(this_indx, :) = opt_aif_c(this_indx, :); %restore the curve to the list of optimal ones
                opt_mean_aif_c = nanmean(opt_aif_c_copy, 1); %update the mean
                set(opt_mean_aif_plot, 'XData', t, 'YData', opt_mean_aif_c) %and update the mean plot
                drawnow     %now                
            end
        else %if clicked pixel is not one of our optimal ones
        end % do nothing
   
        mean_aif_c = opt_mean_aif_c;
        aif_area = trapz(t, mean_aif_c);
    end
    
    else
         mean_aif_c = opt_mean_aif_c;
        aif_area = trapz(t, mean_aif_c);
        
        if key_pushed == 114 || key_pushed == 32 %r is 114 and spacebar is 32; if any of these are pressed, prepare for ROI redrawing
            delete(h_roi)
            for plt = 1:numaifs
                if ishandle(opt_roi_plots(plt)); delete(opt_roi_plots(plt)); end
                if ishandle(opt_aif_plots(plt)); delete(opt_aif_plots(plt)); end
            end
            hold(handles.axes1, 'off')
            hold(handles.axes2, 'off')
%                 wrap_text = ['Click on image to begin drawing of freehand ROI for AIF selection.  '...  
%                 'Program will automatically suggest AIFs in selected ROI.   Use slider to change displayed slice.' ];
%                 set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
%                 blink_this(handles.edit7, 'b') 
        end
        if key_pushed == 115 || key_pushed == 13 %115 is s and 13 is enter
                [file,path] = uiputfile(file_filter);
                if file == 0
                    wrap_text = 'AIF not saved. Selection is still active.';
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                    pause(2)
                else
                path_to_aif = fullfile(path,file);
                try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
                save_data.this_folder = path_to_aif; save_data.this_format = 0;
                save_data.this_name = 0; save_data.data_to_save = mean_aif_c; save_this_file(save_data)
                saving_now = true;
                end
        end
        if key_pushed == 27
            mean_aif_c = false;
            aif_area = false;
            return
        end
        if key_pushed == 99 %c
           mean_aif_c = opt_mean_aif_c; aif_area = trapz(t, mean_aif_c);
            if ishandle(h_roi); delete(h_roi); end
            for plt = 1:numaifs
                if ishandle(opt_roi_plots(plt)); delete(opt_roi_plots(plt)); end
                if ishandle(opt_aif_plots(plt)); delete(opt_aif_plots(plt)); end    
            end
            hold(handles.axes1, 'off')
            hold(handles.axes2, 'off')
            [file,path] = uiputfile(file_filter);
             if file == 0
                    wrap_text = 'AIF not saved. Analysis proceeding...';
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                    pause(3)
            else
            path_to_aif = fullfile(path,file);
            try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
            save_data.this_folder = path_to_aif; save_data.this_format = 0;
            save_data.this_name = 0; save_data.data_to_save = mean_aif_c; save_this_file(save_data)
             end
             return
        end
    end
    end
%     catch %if any error occurs, do not terminate AIF selection
%     end
        %prepare for redrawing on next iteration
        try
if ishandle(h_roi); delete(h_roi); end
if ~saving_now; if ishandle(opt_mean_aif_plot); delete(opt_mean_aif_plot); end; end
for plt = 1:numaifs
    if ishandle(opt_roi_plots(plt)); delete(opt_roi_plots(plt)); end
    if ishandle(opt_aif_plots(plt)); delete(opt_aif_plots(plt)); end    
end
hold(handles.axes1, 'off')
hold(handles.axes2, 'off')
wrap_text = ['Click on image to begin drawing of freehand ROI for AIF selection.   Program will automatically suggest AIFs in selected ROI.   ' ...
    ' Press C to continue with currently displayed AIF.  Press Esc to terminate AIF selection. '];
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
blink_this(handles.edit7, 'b')
        catch
        end
end

end