function cbv_without_mtt = cbv_no_mtt(data, config)

cprintf('*blue', 'CALCULATING CBV: (CBV = k*area{C(t)}/area{AIF}) \n')
cprintf('*blue', '\n')

te = config.te;
tr = config.tr;
img_size = config.img_size;
kappa = config.kappa;
baseline_idx = config.baseline_index;
dsc_data = data.dsc_data;
mask = data.mask;
mean_aif_c = data.mean_aif_c;
save_results = config.save_results;
slice_range = config.slice_range;
target_folder = config.target_folder;

t = 0:tr:(img_size(4)-1)*tr; 

cbv_without_mtt = zeros(img_size(1:3)); %cerebral blood volume

parfor_progress(img_size(1));
for x = slice_range(1):slice_range(2)%1:img_size(1)    
    for y = slice_range(3):slice_range(4)%1:img_size(2)  
        for z = slice_range(5):slice_range(6)%1:img_size_3 
            if mask(x,y,z)
                tmp_tissue_signal   = squeeze(dsc_data(x,y,z,:));
                tmp_tissue_s0       = mean(tmp_tissue_signal(baseline_idx));
                ydata               = -(1/(te)) .* log(tmp_tissue_signal./tmp_tissue_s0);
                % cbv and mtt (method 1)
                cbv_without_mtt(x,y,z) = kappa * trapz(t,ydata)/trapz(t,mean_aif_c).*100; %100 to convert to ml/100g                                              
            end
        end
    end
    parfor_progress;
end
parfor_progress(0);


if save_results

niftiwrite(cbv_without_mtt, strcat(target_folder, '/cbv_without_mtt.nii'))

end

cprintf('magenta','Done. \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')

end