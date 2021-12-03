function display_delay(params_row, row_num, slice_num,caller, handles)
%        This function accepts a row of Delay data, row number, slice number,
%        caller (which deconvolution algorithm is running) and the handles
%        structure. It displays Delay on the GUI axes during deconvolution.
%        This function has no outputs.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slice_range = handles.slice_range;
delay_volume = zeros([128 128 20]);
if strcmp(caller, 'bzd')
    if ~handles.with_delay
        return 
    else
    delay_row = squeeze(params_row(1,:,7)); 
    delay_volume(row_num, :, slice_num) = delay_row;
    image_data = delay_volume;
        current_image = getimage(handles.delay_axes);
        set(handles.delay_image, 'CData', (current_image + squeeze(image_data(:,:,slice_num))))
        drawnow
        if row_num == slice_range(2)
            set(handles.delay_image, 'CData', zeros(128,128))
        end
    end
end

if strcmp(caller, 'svd')
    delay_row = squeeze(params_row(:));
    delay_volume(row_num, :, slice_num) = delay_row;
    image_data = delay_volume;
        current_image = getimage(handles.delay_axes);
        set(handles.delay_image, 'CData', ( current_image + squeeze(image_data(:,:,slice_num))) )
        drawnow
        if row_num == slice_range(2)
            set(handles.delay_image, 'CData', zeros(128,128))
        end
end