function [cbv_BzD_mtt, cbv_oSVD_mtt] = cbv_using_mtt_gui(handles)
%       This function accepts the handles structure and returns CBV
%       calculated as CBF*MTT, for Bezier curve deconvolution and SVD
%       deconvolution.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%display this on main window
% wrap_text = 'Calculating CBV as: CBV = CBF*MTT...';
% set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
%retrieve needed data______________________________________________________
img_size = handles.img_size;
mask = handles.mask;
fitd_cbf = handles.fitd_cbf;
fitd_cbf_svd = handles.fitd_cbf_svd;
mtt_BzD_no_cbv = handles.mtt_BzD_no_cbv;
mtt_oSVD_no_cbv = handles.mtt_oSVD_no_cbv ;
BzD = handles.BzD;
do_SVD = handles.do_SVD;
save_results = handles.save_results;
slice_range = handles.slice_range;

%initialise arrays_________________________________________________________
if BzD; cbv_BzD_mtt = zeros(img_size(1:3)); else; cbv_BzD_mtt = 0;end
if do_SVD; cbv_oSVD_mtt = zeros(img_size(1:3));else; cbv_oSVD_mtt = 0; end
%__________________________________________________________________________
for x = slice_range(1):slice_range(2)
    for y = slice_range(3):slice_range(4)
        for z = slice_range(5):slice_range(6)
            if mask(x,y,z)
                if BzD
                cbv_BzD_mtt(x,y,z) =  fitd_cbf(x,y,z) * mtt_BzD_no_cbv(x,y,z)./60; %/60 to convert to ml/100g
                end
                if do_SVD
                cbv_oSVD_mtt(x,y,z)   =  fitd_cbf_svd(x,y,z) * mtt_oSVD_no_cbv(x,y,z)./60;   
                end
            end
        end
    end
end

%save results______________________________________________________________
if save_results
if BzD
    try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
    save_data.this_folder = handles.target_folder; save_data.this_format = handles.save_format;
    save_data.this_name = 'cbv_as_cbf_times_mtt_BzD'; save_data.data_to_save = cbv_BzD_mtt; save_this_file(save_data)
end
if do_SVD
    try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
    save_data.this_folder = handles.target_folder; save_data.this_format = handles.save_format;
    save_data.this_name = 'cbv_as_cbf_times_mtt_SVD'; save_data.data_to_save = cbv_oSVD_mtt; save_this_file(save_data)
end
end

end