function [mtt_BzD_no_cbv, mtt_oSVD_no_cbv] = mtt_no_CBV_gui(handles)
%       This function accepts the handles structure and returns MTT
%       calculated as area under the residue function, for Bezier curve 
%       deconvolution and deconvolution.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%display this message on the main window
% wrap_text = 'Calculating MTT as area under the residue function...';
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
if BzD; mtt_BzD_no_cbv = zeros(img_size(1:3)); else; mtt_BzD_no_cbv = 0; end
if do_SVD; mtt_oSVD_no_cbv = zeros(img_size(1:3)); else; mtt_oSVD_no_cbv = 0; end
%__________________________________________________________________________
for x = slice_range(1):slice_range(2)
    for y = slice_range(3):slice_range(4)
        for z = slice_range(5):slice_range(6)
            if mask(x,y,z) 
                if BzD
                r_BzD = bezier_residue_function(fitd_omega(x,y,z,:),tv);
                mtt_BzD_no_cbv(x,y,z)   = trapz(tv,r_BzD); 
                end
                if do_SVD
                tmp_r_svd       = squeeze(fitd_r_svd(x,y,z,:)); 
                mtt_oSVD_no_cbv(x,y,z)   = trapz(t,tmp_r_svd)./max(tmp_r_svd(:));
                end
            end
        end
    end
end
%save results______________________________________________________________
if save_results
try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
save_data.this_folder = target_folder; save_data.this_format = handles.save_format;
if BzD; save_data.this_name = 'mtt_as_area_under_r_BzD'; save_data.data_to_save = mtt_BzD_no_cbv; save_this_file(save_data); end
if do_SVD; save_data.this_name = 'mtt_as_area_under_r_SVD'; save_data.data_to_save = mtt_oSVD_no_cbv; save_this_file(save_data); end
end

end