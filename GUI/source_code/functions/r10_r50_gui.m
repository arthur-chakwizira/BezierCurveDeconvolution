function [r10_BzD, r50_BzD, r10_svd, r50_svd] = r10_r50_gui(handles)
%       This function accepts the handles structure and returns R10 and R50
%       which correspond to time points at which the residue function has 
%       fallen to 10% and 50% of its maximum value, respectively.
%       Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%display this on main window
% wrap_text = 'Calculating R10 and/or R50...';
% set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
%retrieve needed data______________________________________________________
tr = handles.tr;
img_size = handles.img_size;
mask = handles.mask;
fitd_omega = handles.fitd_omega;
fitd_r_svd = handles.fitd_r_svd;
BzD = handles.BzD;
do_SVD = handles.do_SVD;
save_results = handles.save_results;
slice_range = handles.slice_range;
target_folder = handles.target_folder;
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
            if mask(x,y,z)
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
end

%save results______________________________________________________________
if save_results
try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
save_data.this_folder = target_folder; save_data.this_format = handles.save_format;
if handles.BzD; save_data.this_name = 'r10_BzD'; save_data.data_to_save = r10_BzD; save_this_file(save_data); end
if handles.do_SVD; save_data.this_name = 'r10_SVD'; save_data.data_to_save = r10_svd; save_this_file(save_data); end
if handles.BzD; save_data.this_name = 'r50_BzD'; save_data.data_to_save = r50_BzD; save_this_file(save_data); end
if handles.do_SVD; save_data.this_name = 'r50_SVD'; save_data.data_to_save = r50_svd; save_this_file(save_data); end
end

end