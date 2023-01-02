function [r_svd, cbf_svd, delay_svd] = SVD_deconvolution_gui(handles)
%        This function handles SVD deconvolution. It does not display
%        results during computation. Its input is the handles structure and
%        its output is: residue functions, cbf and delay
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%retrieve needed data______________________________________________________
te = handles.te;
tr = handles.tr;
img_size = handles.img_size;
save_fitted_params = handles.save_results;
target_folder = handles.target_folder;
slice_range = handles.slice_range;
kappa = handles.kappa;

do_oSVD = handles.do_oSVD;
if do_oSVD; OIndex = handles.OIndex; else; OIndex = 0.095; end

do_cSVD = handles.do_cSVD;
do_sSVD = handles.do_sSVD;
if do_sSVD || do_cSVD; Psvd = handles.Psvd; else; Psvd = 0.2; end

dsc_data = handles.dsc_data;
mask = handles.mask;
mean_aif_c = handles.mean_aif_c;

img_size_23 = img_size(2:3);
t = 0:tr:(img_size(4)-1)*tr;
tl = length(t); 

tmp_baseline_idx = handles.baseline_index;
%__________________________________________________________________________

%initialise arrays_________________________________________________________
r_svd = zeros(img_size); 
cbf_svd = zeros(img_size(1:3));
delay_svd = zeros(img_size(1:3));
%__________________________________________________________________________
%range of deconvolution
xrange = slice_range(1):slice_range(2);
yrange = slice_range(3):slice_range(4);
zrange = slice_range(5):slice_range(6);
%initialise progress bar
% parfor_progress_gui(handles, 'svd', xrange(end)-xrange(1)); 
tic
%__________________________________________________________________________
parfor x = xrange
    tmp_r = zeros([img_size_23 tl]);
    tmp_cbf = zeros(img_size_23);
    tmp_delay = zeros(img_size_23);
    
    oscFact = OIndex;
    Psvd_start = 0.01;
    tmp_t = t; 
    
    for y = yrange
        for z = zrange
            if mask(x,y,z)
                tmp_tissue_signal   = squeeze(dsc_data(x,y,z,:));
                tmp_tissue_s0       = mean(tmp_tissue_signal(tmp_baseline_idx)); 
                C               = -(1/(te)) .* log(tmp_tissue_signal./tmp_tissue_s0);
                    if do_oSVD; [Rut, ~, ~] = oSVD_residue_function(C,mean_aif_c,tmp_t,Psvd_start,oscFact); end
                    if do_cSVD; Rut = cSVD_residue_function(C, mean_aif_c, tmp_t, Psvd); end
                    if do_sSVD; Rut = sSVD_residue_function(C, mean_aif_c, tmp_t, Psvd); end
                    
                     [max_r, tmp_max_id]    = max(Rut);                            
                     tmp_r(y,z,:)    = Rut./max_r;
                    tmp_cbf(y,z) = max_r;%kappa*max(Rut(:))*trapz(tmp_t, C)/( trapz(tmp_t, Rut)*trapz(tmp_t, mean_aif_c));
                    tmp_delay(y,z) = tmp_t(tmp_max_id); 
            end
        end
    end
        r_svd(x,:,:,:) = tmp_r; 
        cbf_svd(x,:,:) = tmp_cbf.*6000; %ml/100g/min
        delay_svd(x,:,:) = tmp_delay; 
%     parfor_progress_gui(handles, 'svd');
end
% parfor_progress_gui(handles, 'svd',0); %close progress bar
 toc
%save data_________________________________________________________________
if save_fitted_params
try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
save_data.this_folder = target_folder; save_data.this_format = handles.save_format;
save_data.this_name = 'residue_functions_SVD'; save_data.data_to_save = r_svd; save_this_file(save_data); 
save_data.this_name = 'delay_SVD'; save_data.data_to_save = delay_svd; save_this_file(save_data); 
end
%__________________________________________________________________________

%report completion
if handles.do_oSVD; wrap_text = 'oSVD deconvolution complete'; end
if handles.do_cSVD; wrap_text = 'cSVD deconvolution complete'; end
if handles.do_sSVD; wrap_text = 'sSVD deconvolution complete'; end
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
end