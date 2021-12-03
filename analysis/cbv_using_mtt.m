function [cbv_BzD_mtt, cbv_oSVD_mtt] = cbv_using_mtt(data,  config)

cprintf('*blue', 'CALCULATING CBV: (CBV = CBF*MTT) \n')
cprintf('*blue', '\n')

img_size = config.img_size;
mask = data.mask;
fitd_cbf = data.fitd_cbf;
fitd_cbf_svd = data.fitd_cbf_svd;
mtt_BzD_no_cbv = data.mtt_BzD_no_cbv;
mtt_oSVD_no_cbv = data.mtt_oSVD_no_cbv ;
BzD = config.BzD;
oSVD = config.oSVD;
save_results = config.save_results;
slice_range = config.slice_range;
target_folder = config.target_folder;



cbv_BzD_mtt = zeros(img_size(1:3)); %cerebral blood volume
cbv_oSVD_mtt = zeros(img_size(1:3));

parfor_progress(img_size(1));
for x = slice_range(1):slice_range(2)%1:img_size(1)  
    for y = slice_range(3):slice_range(4)%1:img_size(2)  
        for z = slice_range(5):slice_range(6)%1:img_size_3
            if mask(x,y,z)
                if BzD
                cbv_BzD_mtt(x,y,z) =  fitd_cbf(x,y,z) * mtt_BzD_no_cbv(x,y,z)./60; %/60 to convert to ml/100g
                end
                if oSVD
                cbv_oSVD_mtt(x,y,z)   =  fitd_cbf_svd(x,y,z) * mtt_oSVD_no_cbv(x,y,z)./60;   
                end
            end
        end
    end
    parfor_progress;
end
parfor_progress(0);


if save_results

niftiwrite(cbv_BzD_mtt, strcat(target_folder ,'/cbv_BzD_mtt.nii'))
niftiwrite(cbv_oSVD_mtt, strcat(target_folder ,'/cbv_oSVD_mtt.nii'))


end

cprintf('magenta','Done. \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')


end