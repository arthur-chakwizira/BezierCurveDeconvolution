function [ttp_from_data, ttp_fit_BzD, ttp_fit_oSVD] = find_ttp_gui(handles)
%       This function accepts the handles structure and returns TTP
%       calculated from input data, and from fitted tissue concentration
%       curve, for Bezier curve deconvolution and SVD deconvolution.
%       Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%No display message; this function runs without consulting the user.

%retrieve needed data______________________________________________________
te = handles.te;
tr = handles.tr;
img_size = handles.img_size;
kappa = handles.kappa;
baseline_idx = handles.baseline_index;
dsc_data = handles.dsc_data;
mask = handles.mask;
mean_aif_c = handles.mean_aif_c;
fitd_cbf = handles.fitd_cbf;
fitd_cbf_svd = handles.fitd_cbf_svd;
fitd_omega = handles.fitd_omega;
fitd_r_svd = handles.fitd_r_svd;
fitd_dk = handles.fitd_dk;
fitd_delay = handles.fitd_delay;
fitd_delay_svd = handles.fitd_delay_svd;
with_dispersion = handles.with_dispersion;
with_delay = handles.with_delay;
BzD = handles.BzD;
do_SVD = handles.do_SVD;
save_results = handles.save_results;
slice_range = handles.slice_range;
t = 0:tr:(img_size(4)-1)*tr;
dt = t(2)-t(1);

%initialise arrays_________________________________________________________
ttp_from_data = zeros(img_size(1:3));
if BzD; ttp_fit_BzD = zeros(img_size(1:3)); else; ttp_fit_BzD = 0; end
if do_SVD; ttp_fit_oSVD = zeros(img_size(1:3)); else; ttp_fit_oSVD = 0; end
%__________________________________________________________________________
for x = slice_range(1):slice_range(2)
    for y = slice_range(3):slice_range(4)
        for z = slice_range(5):slice_range(6)
            if mask(x,y,z) 
                tmp_tissue_signal   = squeeze(dsc_data(x,y,z,:));
                tmp_tissue_s0       = mean(tmp_tissue_signal(baseline_idx));
                ydata               = -(1/(te)) .* log(tmp_tissue_signal./tmp_tissue_s0);      
                % estimate ttp
                [~,c_ttp]           = max(ydata); 
                ttp_from_data(x,y,z)        = t(c_ttp);  
                if BzD
                 r = bezier_residue_function(fitd_omega(x,y,z,:), t); 
                 if with_delay
                 mean_aif_c_BzD = interp1(t,mean_aif_c,t-fitd_delay(x,y,z,:),'pchip','extrap');
                 else
                     mean_aif_c_BzD = mean_aif_c;
                 end
                if with_dispersion
                    dispersed_aif_BzD = dt*filter(mean_aif_c_BzD,1,squeeze(fitd_dk(x,y,z,:)));
                    ct_fit_BzD         = kappa*fitd_cbf(x,y,z,:)*dt.*filter(dispersed_aif_BzD,1,r);
                else
                    ct_fit_BzD = kappa*fitd_cbf(x,y,z,:)*dt.*filter(mean_aif_c_BzD, 1, r);
                end
                 [~,c_ttp_BzD]       = max(ct_fit_BzD);
                 ttp_fit_BzD(x,y,z)  	= t(c_ttp_BzD); 
                end
                
                if do_SVD
                r_svd = fitd_r_svd(x,y,z,:);
                mean_aif_c_oSVD = interp1(t,mean_aif_c,t-fitd_delay_svd(x,y,z,:),'pchip','extrap');
                ct_fit_oSVD = kappa*fitd_cbf_svd(x,y,z,:).*dt.*filter(mean_aif_c_oSVD, 1, r_svd);
                
                [~,c_ttp_oSVD]       = max(ct_fit_oSVD);
                ttp_fit_oSVD(x,y,z)  	= t(c_ttp_oSVD);
                end
                 
            end
        end
    end
end

%save results______________________________________________________________
if save_results
   try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
   save_data.this_folder = handles.target_folder; save_data.this_format = handles.save_format;
   save_data.this_name = 'ttp_without_deconvolution'; save_data.data_to_save = ttp_from_data; save_this_file(save_data)
   if BzD; save_data.this_name = 'ttp_from_fitted_signal_BzD'; save_data.data_to_save = ttp_fit_BzD; save_this_file(save_data); end
   if do_SVD; save_data.this_name = 'ttp_from_calculated_signal_SVD'; save_data.data_to_save = ttp_fit_oSVD; save_this_file(save_data); end
end

end