function [ttp_from_data, ttp_fit_BzD, ttp_fit_oSVD] = find_ttp(data, config)

cprintf('*blue', 'CALCULATING TTP FROM DSC_DATA, BZD-FITTED C(t) AND oSVD-FITTED C(t) \n')
cprintf('*blue', '\n')

te = config.te;
tr = config.tr;
img_size = config.img_size;
kappa = config.kappa;
baseline_idx = config.baseline_index;
dsc_data = data.dsc_data;
mask = data.mask;
mean_aif_c = data.mean_aif_c;
fitd_cbf = data.fitd_cbf;
fitd_cbf_svd = data.fitd_cbf_svd;
fitd_omega = data.fitd_omega;
fitd_r_svd = data.fitd_r_svd;
fitd_dk = data.fitd_dk;
fitd_delay = data.fitd_delay;
fitd_delay_svd = data.fitd_delay_svd;
with_dispersion = config.with_dispersion;
with_delay = config.with_delay;
BzD = config.BzD;
oSVD = config.oSVD;
save_results = config.save_results;
slice_range = config.slice_range;
target_folder = config.target_folder;

t = 0:tr:(img_size(4)-1)*tr; 
dt = t(2)-t(1); 


ttp_from_data = zeros(img_size(1:3)); %cerebral blood volume
ttp_fit_BzD = zeros(img_size(1:3));
ttp_fit_oSVD = zeros(img_size(1:3));


 parfor_progress(img_size(1));
for x = slice_range(1):slice_range(2)%1:img_size(1)   
    for y = slice_range(3):slice_range(4)%1:img_size(2)  
        for z = slice_range(5):slice_range(6)%1:img_size_3 
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
                
                if oSVD
                r_svd = fitd_r_svd(x,y,z,:);
                mean_aif_c_oSVD = interp1(t,mean_aif_c,t-fitd_delay_svd(x,y,z,:),'pchip','extrap');
                ct_fit_oSVD = kappa*fitd_cbf_svd(x,y,z,:).*dt.*filter(mean_aif_c_oSVD, 1, r_svd);
                
                [~,c_ttp_oSVD]       = max(ct_fit_oSVD);
                ttp_fit_oSVD(x,y,z)  	= t(c_ttp_oSVD);
                end
                 
            end
        end
    end
    parfor_progress;
end
parfor_progress(0);


if save_results

niftiwrite(ttp_from_data, strcat(target_folder, '/ttp_from_data.nii'))
niftiwrite(ttp_fit_BzD, strcat(target_folder, '/ttp_fit_BzD.nii'))
niftiwrite(ttp_fit_oSVD,  strcat(target_folder, '/ttp_fit_oSVD.nii'))


end

cprintf('magenta','Done. \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')


end