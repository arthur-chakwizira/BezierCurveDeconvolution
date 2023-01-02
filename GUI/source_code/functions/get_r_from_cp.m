function r = get_r_from_cp(handles)
%       This function accepts the handles structure and returns a 4D array
%       of residue functions from Bezier curve deconvolution. It converts
%       the control points saved in handles to complete residue functions.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%retrieve needed data______________________________________________________
fitd_omega = handles.fitd_omega;
tr = handles.tr;
img_size = handles.img_size;
t = 0:tr:(img_size(4)-1)*tr;
r = zeros(img_size);
slice_range = handles.slice_range;
mask = handles.mask;
%__________________________________________________________________________
%range of deconvolution
xrange = slice_range(1):slice_range(2);
yrange = slice_range(3):slice_range(4);
zrange = slice_range(5):slice_range(6);
%__________________________________________________________________________
for x = xrange
    for y = yrange
        for z = zrange
            if mask(x,y,z)
                r(x,y,z,:) = bezier_residue_function(fitd_omega(x,y,z,:), t); 
            end
        end
    end
end
end