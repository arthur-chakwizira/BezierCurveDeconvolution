function [mean_dsc_signal, t_min_signal] = whole_brain_curve_gui(handles)
%       This function accepts the handles structure and returns a
%       whole-brain signal curve and the time point at which minimum signal
%       occurs.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dsc_data = handles.dsc_data;
mask = handles.mask;
tr = handles.tr; 
img_size = handles.img_size;

t = 0:tr:(img_size(4)-1)*tr; 

tmp_dsc_data = reshape(dsc_data,[prod(img_size(1:3)) img_size(4)]); 
tmp_mask = repmat(double(mask(:)),[1 img_size(4)]);
tmp_mask(tmp_mask==0) = NaN; 
mean_dsc_signal = nanmean(tmp_dsc_data.*tmp_mask,1); 

if isfield(handles, 'plot_whole_brain_curve')
    make_plots = handles.plot_whole_brain_curve;
else
    make_plots =true;
end
    

if make_plots
    set(handles.slider1, 'Visible', 'off')
plot(t,mean_dsc_signal,'k', 'Parent', handles.axes1); xlabel(handles.axes1, 't [s]'); 
ylabel(handles.axes1, 'Signal intensity')
title(handles.axes1, 'Whole brain signal curve')
end
[~,t_min_signal] = min(mean_dsc_signal);
end
