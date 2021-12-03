function [dsc_data, mask] = load_dsc_data(path_to_data, config)


make_mask = config.make_mask;
save_mask = config.save_mask;
mask_th = config.mask_threshold; 

try %if user has specified path to data, use that path
    dsc_nii_file_path = path_to_data;
    dsc_data = niftiread(dsc_nii_file_path);
catch   %otherwise send termination signal
    cprintf('*red','No valid path specified. Analysis will terminate. \n')
    cprintf('*red','\n')
    dsc_data = false; mask = false;
    return
end


% load the data
%niftiread reads files in the nifti format (Neuroimaging Informatics Technology Initiative)
%returns volumetric data in input file
%header information can be obtained using niftiinfo(filename)
dsc_data = double(dsc_data);
img_size = size(dsc_data);

% do a simple masking
[file_path, file_name, ~] = fileparts(dsc_nii_file_path);
path_to_mask = fullfile(file_path, [file_name '_mask.nii']);

if make_mask
    if isfile(path_to_mask)
        mask = niftiread(path_to_mask);
        save_mask = false;
    else
        mask_ref = mean(dsc_data,4); %
        mask = false(img_size(1:3)); %
        mask(mask_ref>mask_th) = 1; %
        se = strel('disk',2);  %structuring element, disk-shaped, radius of 2
        mask = imopen(mask,se); %morphological operator; remove scattered pixels
        mask = imclose(mask,se); %morphological operator; close holes
    end
else
    mask = false;
end

if save_mask
    niftiwrite(double(mask), path_to_mask) % save mask in nifti format; header file is created by default
end

end