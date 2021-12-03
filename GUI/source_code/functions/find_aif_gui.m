function [mean_aif_c, aif_area] = find_aif_gui(handles)
%       This function accepts the handles structure and returns a mean AIF
%       and the area under that mean AIF. It displays a DSC-MRI
%       concentration image on the GUI window and lets the user select 
%       individual pixels in search for an AIF.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dsc_data = handles.dsc_data;
te = handles.te; %s
tr =  handles.tr; %s 
img_size = handles.img_size;
baseline_idx = handles.baseline_index;
mask = handles.mask;

aif_handles = handles.aif_handles;
aif_slice = aif_handles.aif_slice;
aif_time_point = aif_handles.aif_time_point;

hObject = handles.hObject;
file_filter = handles.file_filter;
t = 0:tr:(img_size(4)-1)*tr; 
    
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
mean_aif_c = 0;
aif_area = 0;
%-------------------------------------------------------------------------
no_error = false;
for rep = 1:128^2
    if no_error; break; end
      cla(handles.axes2, 'reset')
      set(handles.axes2,'XTick', [], 'YTick', [])
      set(handles.axes2, 'box', 'on')
    try
set(handles.pushbutton14, 'Visible', 'on')

handles.data_is_for_aif = true;
clim_im = [0 100];
handles.aif_conc_image = imshow(zeros(img_size(1), img_size(2)), clim_im, 'colormap', gray, 'Parent', handles.axes1);
handles.aif_conc_data = dsc_data_c;
set(handles.aif_conc_image, 'CData', squeeze(handles.aif_conc_data(:,:,aif_slice,aif_time_point))  )
title(handles.axes1, ['Concentration: Slice number ' num2str(aif_slice)])
set(handles.slider1, 'Min', 1, 'Max', handles.img_size(3), ...
                         'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', aif_slice)
clim_max = double(max(handles.aif_conc_data(isfinite(handles.aif_conc_data)))) + 1E-10;
if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/999; end
if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end 
set(handles.slider6, 'Visible', 'on') %color setting
set(handles.slider8, 'Visible', 'on')%horizontal scrolling on axes1
set(handles.slider8, 'Min', 1, 'Max', handles.img_size(4), ...
                         'SliderStep', [1, 1]/(handles.img_size(4) - 1), 'Value', aif_time_point)
set(handles.slider6, 'Min', 10E-10, 'Max', clim_max, ...
                         'SliderStep', [step_small, step_big], 'Value', clim_im(2))
set(handles.slider10, 'Visible', 'on') %color setting
set(handles.slider10, 'Min', 0, 'Max', clim_max-1E-10, ...
                         'SliderStep', [step_small, step_big], 'Value', 0)                     
guidata(hObject, handles)
wrap_text = 'Use slider to select slice for AIF selection. Push Continue when done.';
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
blink_this(handles.edit7, 'b')
uiwait;

handles.data_is_for_aif = false;
set(handles.pushbutton14, 'Visible', 'off')
guidata(hObject, handles)
hold(handles.axes1, 'on')

wrap_text = ['Click on image to select AIFs. Click once to activate cursor.'...
    'Clicking on selected pixel will deselect it.   Press:  S/Enter to save, R/Spacebar to redo selection, C to continue with current AIF,'...
    ' Esc to terminate selection.'];
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
blink_this(handles.edit7, 'b')    
    
    n_clicks = 128*128; %maximum number of clicks is this. We don't expect more :-)
    current_slice = round(get(handles.slider1, 'Value'));
    handles.coordinates = zeros(n_clicks,2); %xi, yi
    handles.roi_plots = zeros(n_clicks, 1);
    handles.aif_c = NaN(n_clicks, img_size(4));
    handles.aif_plots = zeros(n_clicks,1);
    handles.xcoords = zeros(n_clicks,1);
    handles.ycoords = zeros(n_clicks,1);
    handles.mean_aif_plot = plot(0, 0, 'k-','LineWidth', 2, 'Parent', handles.axes2);
%     xlim(handles.axes2, [0 img_size(4)])
%     ylim(handles.axes2, [-10 110])
%     set(handles.axes2, 'YTick', 0:10:100)
    xlabel(handles.axes2, 'time [s]')
    ylabel(handles.axes2, 'Concentration')
    title(handles.axes2, 'AIF selection')
    legend(handles.axes2, 'Mean')
    hold(handles.axes2, 'on')
    for i = 1:n_clicks
        exit_now = waitforbuttonpress; %0 for mouse-click, 1 for key press
        clear('this_key')
        this_key = double(get(handles.figure1, 'CurrentCharacter'));
        if exit_now
            if this_key == 115 || this_key == 13 %s or enter
                mean_aif_c = handles.mean_aif_c;
                aif_area = trapz(t, mean_aif_c);
                [file,path] = uiputfile(file_filter);
                if file ~= 0
                path_to_aif = fullfile(path,file);
                try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
                save_data.this_folder = path_to_aif; save_data.this_format = 0;
                save_data.this_name = 0; save_data.data_to_save = mean_aif_c; save_this_file(save_data)       
                wrap_text = ['Mean AIF saved. AIF selection still active. '...
                    'Press:  S/Enter to save, R/Spacebar to redo selection, C to continue with current AIF,'...
                    ' Esc to terminate selection.'];
                set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
                blink_this(handles.edit7, 'b')
                else
                    wrap_text = ['Mean AIF not saved. AIF selection still active. '...
                    'Press:  S/Enter to save, R/Spacebar to redo selection, C to continue with current AIF,'...
                    'Esc to terminate selection.'];
                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
%                     blink_this(handles.edit7, 'r')                    
                end
            end
            if this_key == 114 || this_key == 32 %r or spacebar
                no_error = false;
               break
            end
            if this_key == 27
                mean_aif_c = false; aif_area = false;
                return
            end %escape
            if this_key == 99
                mean_aif_c = handles.mean_aif_c;
                aif_area = trapz(t, mean_aif_c);
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
                wrap_text = 'AIF selection complete.';
                set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
                end
                return
            end %c             
        end
        [x, y] = the_ginput(1);
        pixel_x = round(y); 
        pixel_y = round(x);
        [already_clicked, this_pixel] = check_existence(i);
        if already_clicked && this_pixel ~= false
            set(handles.aif_plots(this_pixel), 'XData', [], 'YData', [])
            set(handles.roi_plots(this_pixel), 'XData', [], 'YData', [])
             handles.aif_c(this_pixel, :) = NaN(1, img_size(4));
             handles.xcoords(this_pixel) = false;
             handles.ycoords(this_pixel) = false;
             drawnow
        else
        handles.xcoords(i) = pixel_x;
        handles.ycoords(i) = pixel_y;
        handles.coordinates(i, :) = [x, y];
               
        handles.aif_c(i,:) = dsc_data_c(pixel_x, pixel_y, current_slice, :);
        handles.aif_plots(i) = plot(t, handles.aif_c(i, :),'LineWidth', 0.8, 'HandleVisibility', 'off', 'Parent', handles.axes2);
        this_color = get(handles.aif_plots(i), 'color');
        handles.roi_plots(i) = plot(handles.coordinates(i,1), handles.coordinates(i,2), 's',...
            'MarkerFaceColor', this_color, 'MarkerEdgeColor', this_color, 'MarkerSize', 4);
        
        drawnow
        end
        handles.mean_aif_c = nanmean(handles.aif_c(1:i, :), 1);
        set(handles.mean_aif_plot, 'XData', t, 'YData', handles.mean_aif_c)
        drawnow
        guidata(hObject, handles)
        no_error = true;
    end
    hold(handles.axes2, 'off')

    catch
        no_error = false;
        cla(handles.axes2, 'reset')
        set(handles.axes2,'XTick', [], 'YTick', [])
        set(handles.axes2, 'box', 'on')
    end

%clean up
handles.coordinates = []; %xi, yi
handles.roi_plots = [];
handles.aif_c = [];
handles.aif_plots = [];
handles.xcoords = [];
handles.ycoords = [];
handles.mean_aif_plot = [];
guidata(hObject, handles)
mean_aif_c = handles.mean_aif_c;
aif_area = trapz(t, mean_aif_c);
end

    function [exist_s, this_indx] = check_existence(indx)
        exist_s = false;
        this_indx = false;
        for c = 1:indx %check is points have already been clicked
            if handles.xcoords(c) == pixel_x & handles.ycoords(c) == pixel_y
                exist_s = true;
                this_indx = c;
            end
        end
    end

end
