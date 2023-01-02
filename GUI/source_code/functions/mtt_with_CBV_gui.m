function [mtt_BzD_cbv, mtt_oSVD_cbv] = mtt_with_CBV_gui(handles)
%       This function accepts the handles structure and returns MTT
%       calculated as CBV/CBF, for Bezier curve deconvolution and SVD
%       deconvolution.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%display this message on the main window
% wrap_text = 'Calculating MTT as MTT = CBV/CBF...';
% set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
%get needed data___________________________________________________________
img_size = handles.img_size;
mask = handles.mask;
cbv_without_mtt = handles.cbv_without_mtt;
fitd_cbf = handles.fitd_cbf;
fitd_cbf_svd = handles.fitd_cbf_svd;
BzD = handles.BzD;
do_SVD = handles.do_SVD;
save_results = handles.save_results;
slice_range = handles.slice_range;
target_folder = handles.target_folder;
%__________________________________________________________________________
if BzD; mtt_BzD_cbv = zeros(img_size(1:3)); else; mtt_BzD_cbv = 0; end
if do_SVD; mtt_oSVD_cbv = zeros(img_size(1:3)); else; mtt_oSVD_cbv = 0; end
%__________________________________________________________________________
for x = slice_range(1):slice_range(2)
    for y = slice_range(3):slice_range(4)
        for z = slice_range(5):slice_range(6)
            if mask(x,y,z)
                if BzD
                    mtt_BzD_cbv(x,y,z) = cbv_without_mtt(x,y,z)/fitd_cbf(x,y,z).*60; %x60 to convert to seconds
                end
                if do_SVD
                    mtt_oSVD_cbv(x,y,z) = cbv_without_mtt(x,y,z)/fitd_cbf_svd(x,y,z).*60; %x60 to convert to seconds
                end
            end
        end
    end
end


%save results______________________________________________________________
if save_results
try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
save_data.this_folder = target_folder; save_data.this_format = handles.save_format;
if BzD; save_data.this_name = 'mtt_as_cbv_over_cbf_BzD'; save_data.data_to_save = mtt_BzD_cbv; save_this_file(save_data); end
if do_SVD; save_data.this_name = 'mtt_as_cbv_over_cbf_SVD'; save_data.data_to_save = mtt_oSVD_cbv; save_this_file(save_data); end
end

end