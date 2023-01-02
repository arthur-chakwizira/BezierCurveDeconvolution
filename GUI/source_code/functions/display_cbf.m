function display_cbf(params_row, row_num, slice_num,caller, handles)
%        This function accepts a row of CBF data, row number, slice number,
%        caller (which deconvolution algorithm is running) and the handles
%        structure. It displays CBF on the GUI axes during deconvolution.
%        This function has no outputs.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slice_range = handles.slice_range;
cbf_volume = zeros([128 128 20]);
if strcmp(caller, 'bzd')
    cbf_row = squeeze(params_row(1,:,6))*6000; 
    cbf_volume(row_num, :, slice_num) = cbf_row;
    image_data = cbf_volume;
        current_image = getimage(handles.cbf_axes);
        set(handles.cbf_image, 'CData', (current_image + squeeze(image_data(:,:,slice_num))))
        drawnow
        if row_num == slice_range(2)
            set(handles.cbf_image, 'CData', zeros(128,128))
        end
end

if strcmp(caller, 'svd')
    cbf_row = squeeze(params_row(:))*6000;
    cbf_volume(row_num, :, slice_num) = cbf_row;
    image_data = cbf_volume;
        current_image = getimage(handles.cbf_axes);
        set(handles.cbf_image, 'CData', ( current_image + squeeze(image_data(:,:,slice_num))) )
        drawnow
        if row_num == slice_range(2)
            set(handles.cbf_image, 'CData', zeros(128,128))
        end
end
