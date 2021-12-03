function [mtt_BzD_cbv, mtt_oSVD_cbv] = mtt_with_CBV(data, config)

cprintf('*blue', 'CALCULATING MTT: (MTT = CBV/CBF) \n')
cprintf('*blue', '\n')

img_size = config.img_size;
mask = data.mask;
cbv_without_mtt = data.cbv_without_mtt;
fitd_cbf = data.fitd_cbf;
fitd_cbf_svd = data.fitd_cbf_svd;
BzD = config.BzD;
oSVD = config.oSVD;
save_results = config.save_results;
slice_range = config.slice_range;
target_folder = config.target_folder;

mtt_BzD_cbv = zeros(img_size(1:3)); %cerebral blood volume
mtt_oSVD_cbv = zeros(img_size(1:3));

 parfor_progress(img_size(1));
for x = slice_range(1):slice_range(2)%1:img_size(1)    
    for y = slice_range(3):slice_range(4)%1:img_size(2) 
        for z = slice_range(5):slice_range(6)%1:img_size_3 
            if mask(x,y,z) 
                if BzD
                    mtt_BzD_cbv(x,y,z) = 60*squeeze(cbv_without_mtt(x,y,z))./squeeze(fitd_cbf(x,y,z)); %x60 to convert to seconds 
                end
                if oSVD
                    mtt_oSVD_cbv(x,y,z) = 60*squeeze(cbv_without_mtt(x,y,z))./squeeze(fitd_cbf_svd(x,y,z)); %x60 to convert to seconds 
                end
            end
        end
    end
    parfor_progress;
end
parfor_progress(0);


if save_results
niftiwrite(mtt_BzD_cbv, strcat(target_folder, '/mtt_BzD_cbv.nii'))
niftiwrite(mtt_oSVD_cbv, strcat(target_folder, '/mtt_oSVD_cbv.nii'))
end

cprintf('magenta','Done. \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')


end