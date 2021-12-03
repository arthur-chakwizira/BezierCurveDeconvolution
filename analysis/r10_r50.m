function [r10_BzD, r50_BzD, r10_svd, r50_svd] = r10_r50(data,config)
cprintf('*blue', 'CALCULATING R10 and R50 \n')
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

r10_BzD = zeros(img_size(1:3)); %cerebral blood volume
r50_BzD = zeros(img_size(1:3));

r10_svd = zeros(img_size(1:3)); %cerebral blood volume
r50_svd = zeros(img_size(1:3));


parfor_progress(img_size(1));
%parfor_progress(11);
for x = slice_range(1):slice_range(2)%1:img_size(1)   
    for y = slice_range(3):slice_range(4)%1:img_size(2)     
        for z = slice_range(5):slice_range(6)%1:img_size_3 
            if mask(x,y,z)
                if BzD
                r = bezier_residue_function(fitd_omega(x,y,z,:), tv);
                
                [~,r10_idx]     = min(abs(r-0.1)); 
                r10_BzD(x,y,z)    = tv(r10_idx);
                
                [~,r50_idx]     = min(abs(r-0.5)); 
                r50_BzD(x,y,z)    = tv(r50_idx); 
                end
                
                if oSVD
                try
                    r_osvd = interp1(t,fitd_r_svd(x,y,z,:),tv,'pchip');
                 [~,r10_idx_svd]     = min(abs(r_osvd-0.1)); 
                 r10_svd(x,y,z)    = tv(r10_idx_svd); 
                     
                 [~,r50_idx_svd]     = min(abs(r_osvd-0.5)); 
                 r50_svd(x,y,z)    = tv(r50_idx_svd); 

                catch
                   r_osvd = fitd_r_svd(x,y,z,:);
                 [~,r10_idx_svd]     = min(abs(r_osvd-0.1));
                 r10_svd(x,y,z)    = t(r10_idx_svd); 
                     
                 [~,r50_idx_svd]     = min(abs(r_osvd-0.5));
                 r50_svd(x,y,z)    = t(r50_idx_svd); 
                end
                end
            end
        end
    end
    parfor_progress;
end
parfor_progress(0);


if save_results
niftiwrite(r10_BzD, strcat(target_folder, '/r10_BzD.nii'))
niftiwrite(r50_BzD, strcat(target_folder, '/r50_BzD.nii'))

end
cprintf('magenta','Done. \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')


end