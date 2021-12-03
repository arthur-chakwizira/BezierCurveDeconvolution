function [mtt_BzD_no_cbv, mtt_oSVD_no_cbv] = mtt_no_CBV(data, config)

cprintf('*blue', 'CALCULATING MTT: (MTT = area{R(t)}) \n')
cprintf('*blue', '\n')

tr = config.tr;
img_size = config.img_size;
mask = data.mask;
fitd_omega = data.fitd_omega;
fitd_r_svd = data.fitd_r_svd;
BzD = config.BzD;
oSVD = config.oSVD;
save_results = config.save_results;
slice_range = config.slice_range;
target_folder = config.target_folder;

t = 0:tr:(img_size(4)-1)*tr; 
tv = 0:0.1:t(end);


mtt_BzD_no_cbv = zeros(img_size(1:3)); %cerebral blood volume
mtt_oSVD_no_cbv = zeros(img_size(1:3));


 parfor_progress(img_size(1));
for x = slice_range(1):slice_range(2)%1:img_size(1)  
    for y = slice_range(3):slice_range(4)%1:img_size(2)  
        for z = slice_range(5):slice_range(6)%1:img_size_3
            if mask(x,y,z) 
                if BzD
                r_BzD = bezier_residue_function(fitd_omega(x,y,z,:),tv);
                mtt_BzD_no_cbv(x,y,z)   = trapz(tv,r_BzD); 
                end
                if oSVD
                tmp_r_svd       = squeeze(fitd_r_svd(x,y,z,:)); 
                mtt_oSVD_no_cbv(x,y,z)   = trapz(t,tmp_r_svd);
                end
            end
        end
    end
    parfor_progress;
end
parfor_progress(0);


if save_results

niftiwrite(mtt_BzD_no_cbv, strcat(target_folder, '/mtt_BzD_no_cbv.nii'))
niftiwrite(mtt_oSVD_no_cbv, strcat(target_folder, '/mtt_oSVD_no_cbv.nii'))

end
cprintf('magenta','Done. \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')


end