function r10_r50_gui_indep(handles)
%       This function accepts the handles structure and returns R10 and R50
%       which correspond to time points at which the residue function has
%       fallen to 10% and 50% of its maximum value, respectively.
%It does the calculation independent of deconvolution and thus
%requires residue functions to be loaded
%       Author:
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wrap_text = 'Select input data location';
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
file_folder = uigetdir('', 'Select input data location.'); %invoke file-select dialog with this title
if file_folder == 0 %if user cancels dialog
    wrap_text = 'No input path selected. Analysis terminated.'; %report and terminate execution
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
    return
end

[base_folder, ~] = fileparts(file_folder);
handles.file_folder = base_folder;
[~, tr] = get_te_tr(handles);


%display this on main window
% wrap_text = 'Calculating R10 and/or R50...';
% set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
%retrieve needed data______________________________________________________
fitd_r_svd = handles.fitd_r_svd;
BzD = handles.BzD;
if BzD
    omega_fn = fullfile(file_folder, 'control_points_BzD.nii');
    if ~isfile(omega_fn);  wrap_text = 'Selected folder does not contain Bezier control points.'; set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r'); return; end
    fitd_omega = niftiread(omega_fn);
    img_size = size(fitd_omega, [1 2 3]);
    slice_range = [1 img_size(1)  1  img_size(2)   1    img_size(3)];
end
do_SVD = handles.do_SVD;
if do_SVD
    r_svd_fn = fullfile(file_folder, 'fitd_r_svd.nii');
    if ~isfile(omega_fn);  wrap_text = 'Selected folder does not contain SVD residue functions.'; set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r'); return; end
    fitd_r_svd = niftiread(r_svd_fn);
    img_size = size(fitd_r_svd, [1 2 3]);
    slice_range = [1 img_size(1)  1  img_size(2)   1    img_size(3)];
end

target_folder = file_folder;
t = 0:tr:(img_size(4)-1)*tr;
tv = 0:0.1:t(end);
%initialise arrays_________________________________________________________
if BzD; r10_BzD = zeros(img_size(1:3)); r50_BzD = zeros(img_size(1:3));
else; r10_BzD = 0; r50_BzD = 0; end
if do_SVD; r10_svd = zeros(img_size(1:3)); r50_svd = zeros(img_size(1:3));
else; r10_svd = 0; r50_svd = 0; end
%__________________________________________________________________________
for x = slice_range(1):slice_range(2)
    for y = slice_range(3):slice_range(4)
        for z = slice_range(5):slice_range(6)
                if BzD
                    r = bezier_residue_function(fitd_omega(x,y,z,:), tv);
                    [~,r10_idx]     = min(abs(r-0.1));
                    r10_BzD(x,y,z)    = tv(r10_idx);
                    
                    [~,r50_idx]     = min(abs(r-0.5));
                    r50_BzD(x,y,z)    = tv(r50_idx);
                end
                
                if do_SVD
                    try
                        r_osvd = interp1(t,fitd_r_svd(x,y,z,:),tv,'pchip');
                        [~,r10_idx_svd]     = min(abs(r_osvd-0.1));
                        r10_svd(x,y,z)    = tv(r10_idx_svd);
                        
                        [~,r50_idx_svd]     = min(abs(r_osvd-0.5));
                        r50_svd(x,y,z)    = tv(r50_idx_svd);
                    catch
                        r_osvd = fitd_r_svd(x,y,z,:);
                        [~,r10_idx_svd]     = min(abs(r_osvd-0.1));
                        r10_svd(x,y,z)    = t(r10_idx_svd);
                        
                        [~,r50_idx_svd]     = min(abs(r_osvd-0.5));
                        r50_svd(x,y,z)    = t(r50_idx_svd);
                    end
                end
        end
    end
end

%save results______________________________________________________________
try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
save_data.this_folder = target_folder; save_data.this_format = handles.save_format;
if handles.BzD; save_data.this_name = 'r10_BzD'; save_data.data_to_save = r10_BzD; save_this_file(save_data); end
if handles.do_SVD; save_data.this_name = 'r10_SVD'; save_data.data_to_save = r10_svd; save_this_file(save_data); end
if handles.BzD; save_data.this_name = 'r50_BzD'; save_data.data_to_save = r50_BzD; save_this_file(save_data); end
if handles.do_SVD; save_data.this_name = 'r50_SVD'; save_data.data_to_save = r50_svd; save_this_file(save_data); end

end