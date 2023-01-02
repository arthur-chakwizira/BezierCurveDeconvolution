function [dsc_data, mask, header_info] = load_dsc_data_gui(path_to_data, handles)
%       This function accepts a path to input data and the handles structure
%       and returns a 4D array of DSC-MRI data, a mask if relevant and the
%       header information of the input file.
%        Author: 
%              Arthur Chakwizira
%              arthur.chakwizira@med.lu.se
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if ~isfile(path_to_data) && ~isfolder(path_to_data)
    wrap_text = 'No valid path specified. Analysis terminated.';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
     opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
     errordlg('\fontsize{10}No valid path specified..', 'Loading DSC-MRI data', opts);
    dsc_data = false; mask = false; header_info = [];
    return
end

dsc_pre_data = read_this_file(path_to_data, 4);
if isstruct(dsc_pre_data) %happens when input file is dicom
    dsc_data = dsc_pre_data.the_data;
    header_info = dsc_pre_data.the_info;
else
    if isequal(dsc_pre_data,false) %happens when data loading failed
    wrap_text = ['Failed to load data from: ' path_to_data];
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
    dsc_data = false; mask = false; header_info = [];
    return
    else
        dsc_data = dsc_pre_data; %happens with everything else
        header_info = [];
    end
end

dsc_data = double(dsc_data);
img_size = size(dsc_data);

% do a simple masking
[file_path, file_name, ~] = fileparts(path_to_data);
path_to_mask = fullfile(file_path, [file_name '_mask.nii']);

mask_th = handles.mask_threshold;

if isfile(path_to_mask)
    mask = niftiread(path_to_mask);
else
    mask_ref = mean(dsc_data,4);
     mask = false(img_size(1:3));
    mask(mask_ref>mask_th) = 1; 
    se = strel('disk',2); 
    mask = imopen(mask,se); 
    mask = imclose(mask,se); 
end


   try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
   save_data.this_folder = file_path; save_data.this_format = handles.save_format;
   save_data.this_name = [file_name '_mask.nii']; save_data.data_to_save = double(mask); save_this_file(save_data)    

end