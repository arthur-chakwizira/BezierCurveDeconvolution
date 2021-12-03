function [r_svd, cbf_svd, delay_svd] = oSVD_deconvolution(data, config)

%Computing residue functions using oSVD 

te = config.te;
tr = config.tr;
img_size = config.img_size;
% notify_when_done = config.notify_when_done;
save_fitted_params = config.save_fitted_params;
target_folder = config.target_folder;
slice_range = config.slice_range;
% kappa = config.kappa;

dsc_data = data.dsc_data;
mask = data.mask;
mean_aif_c = data.mean_aif_c;

img_size_23 = img_size(2:3); 
t = 0:tr:(img_size(4)-1)*tr; 
tl = length(t); 

tmp_baseline_idx = config.baseline_index;

%if no slice range has been specified, calculate all slices
if  isempty(slice_range)
    slice_range = [1 img_size(3)];
end

r_svd = zeros(img_size);
cbf_svd = zeros(img_size(1:3)); 
delay_svd = zeros(img_size(1:3));

xrange = slice_range(1):slice_range(2);
yrange = slice_range(3):slice_range(4);
zrange = slice_range(5):slice_range(6);

parfor_progress(img_size(1));
parfor x = xrange%1:img_size(1) 
           
    tmp_r = zeros([img_size_23 tl]); 
    tmp_cbf = zeros(img_size_23); 
    tmp_delay = zeros(img_size_23);
    
    oscFact = 0.095;%0.065;%0.035;
    Psvd_start = 0.01; %0.1
    tmp_t = t; 
    
    for y = yrange%1:img_size_2
        for z = zrange%1:img_size_3 
            if mask(x,y,z)
                
                tmp_tissue_signal   = squeeze(dsc_data(x,y,z,:));
             
                tmp_tissue_s0       = mean(tmp_tissue_signal(tmp_baseline_idx));  
                C               = -(1/(te)) .* log(tmp_tissue_signal./tmp_tissue_s0);
                       
                    [Rut, ~, ~] = oSVD_residue_function(C,mean_aif_c,tmp_t,Psvd_start,oscFact);
%                     Rut = sSVD(C, mean_aif_c, tmp_t); %Uncomment for standard svd deconvolution  
                     tmp_r(y,z,:)    = Rut; 
                    [tmp_cbf(y,z), tmp_max_id]    = max(Rut);
%                     tmp_cbf(y,z) = kappa*max(Rut(:))*trapz(tmp_t, C)/( trapz(tmp_t, Rut)*trapz(tmp_t, mean_aif_c));
                    tmp_delay(y,z) = tmp_t(tmp_max_id);
            end
        end
    end
        
        r_svd(x,:,:,:) = tmp_r; 
        cbf_svd(x,:,:) = tmp_cbf.*6000; 
        delay_svd(x,:,:) = tmp_delay;
        
    parfor_progress;
end
parfor_progress(0);

    
if save_fitted_params
        niftiwrite(r_svd, strcat(target_folder, '/fitd_r_svd.nii'))
        niftiwrite(cbf_svd, strcat(target_folder, '/fitd_cbf_svd.nii'))
        niftiwrite(delay_svd, strcat(target_folder, '/fitd_delay_svd.nii'))
end

cprintf('magenta','oSVD deconvolution complete. \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')
end