function [mean_dsc_signal, t_min_signal] = whole_brain_curve(data, config)

%plotting signal as a function of time gives quantities such as arrival
%time, negative enhancement integral, mean time to enhance, etc

dsc_data = data.dsc_data;
mask = data.mask;
tr = config.tr; %seconds
img_size = config.img_size;
make_plots = config.make_plots;
x
t = 0:tr:(img_size(4)-1)*tr;


tmp_dsc_data = reshape(dsc_data,[prod(img_size(1:3)) img_size(4)]); %convert dsc_data into 2D matrix
tmp_mask = repmat(double(mask(:)),[1 img_size(4)]); %turn mask into same shape
tmp_mask(tmp_mask==0) = NaN; %convert zeros in mask to undefined
mean_dsc_signal = nanmean(tmp_dsc_data.*tmp_mask,1); %apply mask to dsc_data and compute mean 
%nanmean computes the mean after removing all NaN values. 
%The mean calculated here is the average over all voxels across all slices
%at each time point, resulting in a total of N_time_points data points.
%Hence, "whole brain signal curve"

if make_plots
figure('position',[200 200 1200 800])
plot(t,mean_dsc_signal,'k'); xlabel('t [s]'); ylabel('S(t)') %plot mean whole brain signal as a function of time
xlabel('t [s]'); ylabel('Signal intensity')
title('Whole brain signal curve')
end
[~,t_min_signal] = min(mean_dsc_signal);
end
