function display_mtt(params_row, row_num, slice_num,caller, handles)
%        This function accepts a row of MTT data, row number, slice number,
%        caller (which deconvolution algorithm is running) and the handles
%        structure. It displays MTT on the GUI axes during deconvolution.
%        This function has no outputs.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slice_range = handles.slice_range;
mtt_volume = zeros([128 128 20]);
tr = handles.tr;
img_size = handles.img_size;
t = 0:tr:(img_size(4)-1)*tr;

if strcmp(caller, 'bzd')
    omega_row = squeeze(params_row(1,:,1:5)); 
    bzd_r_row = zeros(size(omega_row, 1), img_size(4)); 
    mtt_row = zeros(1, img_size(2));
    for r = 1:size(bzd_r_row, 1)
        bzd_r_row(r, :) = bezier_residue_function(squeeze(omega_row(r,:)), t);
        mtt_row(r) = trapz(t, squeeze(bzd_r_row(r, :)));
    end
    mtt_volume(row_num, :, slice_num) = mtt_row; 
    image_data = mtt_volume;
        current_image = getimage(handles.mtt_axes);
        set(handles.mtt_image, 'CData', (current_image + squeeze(image_data(:,:,slice_num))))
        drawnow
        if row_num == slice_range(2)
            set(handles.mtt_image, 'CData', zeros(128,128))
        end
end

if strcmp(caller, 'svd')
    svd_r_row = squeeze(params_row(1,:,:)); 
    mtt_row = zeros(1, img_size(2));
    for r = 1:size(svd_r_row, 1)
        r_norm = squeeze(svd_r_row(r, :));
        mtt_row(r) = trapz(t, r_norm)./max(r_norm(:));
    end    
    mtt_volume(row_num, :, slice_num) = mtt_row; 
    image_data = mtt_volume;
        current_image = getimage(handles.mtt_axes);
        set(handles.mtt_image, 'CData', ( current_image + squeeze(image_data(:,:,slice_num))) )
        drawnow
        if row_num == slice_range(2)
            set(handles.mtt_image, 'CData', zeros(128,128))
        end
end
