function cbv_without_mtt = cbv_no_mtt_gui(handles)
%       This function accepts the handles structure and returns CBV
%       calculated as the ratio of the area under the tissue concentration
%       curve to the area under the AIF.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%display this message on the main window
% wrap_text = 'Calculating CBV as area under C(t) divided by area under AIF...)';
% set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
%retrieve needed data______________________________________________________
te = handles.te;
tr = handles.tr;
img_size = handles.img_size;
kappa = handles.kappa;
baseline_idx = handles.baseline_index;
dsc_data = handles.dsc_data;
mask = handles.mask;
mean_aif_c = handles.mean_aif_c;
slice_range = handles.slice_range;
%__________________________________________________________________________
t = 0:tr:(img_size(4)-1)*tr; %time vector
cbv_without_mtt = zeros(img_size(1:3)); %initialise CBV array
%__________________________________________________________________________
for x = slice_range(1):slice_range(2)
    for y = slice_range(3):slice_range(4)
        for z = slice_range(5):slice_range(6)
            if mask(x,y,z) 
                tmp_tissue_signal   = squeeze(dsc_data(x,y,z,:));
                tmp_tissue_s0       = mean(tmp_tissue_signal(baseline_idx)); %signal
                ydata               = -(1/(te)) .* log(tmp_tissue_signal./tmp_tissue_s0); %concentration
                cbv_without_mtt(x,y,z) = kappa * trapz(t,ydata)/trapz(t,mean_aif_c).*100; %100 to convert to ml/100g                                              
            end
        end
    end
end

end