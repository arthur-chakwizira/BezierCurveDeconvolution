function [mean_aif_c, aif_area] = find_aif(data, config)

dsc_data = data.dsc_data;
mask = data.mask;
te = config.te; %s
tr =  config.tr; %s
img_size = config.img_size;
baseline_idx = config.baseline_index;
tail_idx = config.tail_index;

aif_config = config.aif_config;
aif_slice = aif_config.aif_slice;
make_plots = aif_config.make_plots;
save_aif = aif_config.save_aif;
aif_time_point = aif_config.aif_time_point;

t = 0:tr:(img_size(4)-1)*tr; %time points at which each of the DSC-MRI images were taken

% start to search for AIF:s in manually drawn ROI
figure('position',[200 200 1200 800])
subplot(1,2,1)
colormap(gray) 

% show slice
imagesc(squeeze(dsc_data(:,:,aif_slice,aif_time_point))), axis image off 
%This aif_slice is chosen at the time point corresponding to the minimum on the whole-brain signal curve

h_roi = drawfreehand; %create draggable freehand region
aif_roi = h_roi.createMask; %createMask is a function of the imfreehand object. It creates a mask
%within the image: ones in ROI and zeros outside.
aif_mask = false(img_size(1:3));

% pick out curves from the ROI
aif_slices = (aif_slice-1):(aif_slice+1); %include 3 slices
for i = 1:length(aif_slices) 
    aif_mask(:,:,aif_slices(i)) = aif_roi.*mask(:,:,aif_slices(i));
end

tmp_dsc_data = reshape(dsc_data,[prod(img_size(1:3)) img_size(4)]); 
tmp_mask = repmat(aif_mask(:),[1 img_size(4)]); 
aif_data = reshape(tmp_dsc_data(tmp_mask),[sum(aif_mask(:)) img_size(4)]);

% search for aif based on FM and Ymax

if make_plots
    figure(2)
    plot(t,mean(aif_data),'k'), hold on
    xlabel('t [s]'); ylabel('AIF(t)')
    title('mean aif (raw signal)')
end

% convert to concentration
aif_data_c = aif_data.*0; %

% find the indicies with top 5% highest peaks
ymax    = zeros(sum(aif_mask(:)),1);    
fm      = zeros(sum(aif_mask(:)),1);
noise   = zeros(sum(aif_mask(:)),1);
% pointy  = zeros(sum(aif_mask(:)),1);

for i = 1:sum(aif_mask(:))
    
    aif_s0 = mean(mean(aif_data(i,baseline_idx)));

    aif_data_c(i,:) = -(1/te) .* log(aif_data(i,:)./aif_s0) ;
    aif_data(i,:) = aif_data(i,:);
    
    ymax(i) = max(aif_data_c(i,:));
    fm(i) = sum(t.*aif_data_c(i,:))/sum(aif_data_c(i,:));
    
    % noise level
    noise(i) = std(aif_data(i,baseline_idx)) + std(aif_data(i,tail_idx));
    
    % remove to pointy curves
%     pointy(i) = sum(abs(diff(aif_data(i,:))));
end


% high peaks ymax
numaifs = 5;
[~,sortIndex] = sort(ymax(:),'descend');
ymax_index = sortIndex(1:numaifs); % highest values

% remove noisy and pointy
% ymax_index(noise(ymax_index)>40) = [];
% ymax_index(pointy(ymax_index)>2200) = [];

% remove low s0
% ymax_index(aif_s0(ymax_index)<600) = [];


if make_plots
    figure('position',[200 200 1400 400])
    subplot(2,3,1)
    plot(t,aif_data_c(ymax_index,:),'k'), hold on, %ylim([0 120])
    xlabel('t [s]'); ylabel('AIF(t)')
    title('curves with highest values')
end


% low fm
numaifs = round(sum(aif_mask(:))/5);
[~,sortIndex] = sort(fm(:),'ascend');
fm_index = sortIndex(1:numaifs); 

% remove noisy and pointy
%fm_index(noise(fm_index)>10) = [];
%fm_index(pointy(fm_index)>400) = [];

% remove low s0
% ymax_index(aif_s0(ymax_index)<600) = [];

if make_plots
    subplot(2,3,2)
    plot(t,aif_data_c(fm_index,:),'k'); hold on;
    xlabel('t [s]'); ylabel('AIF(t)')
    title('lowest fm values')
end

% find the intersect
try
    ymaxfm_index = intersect(ymax_index,fm_index); %curves with lowest fm and max values
%     numaifs = length(ymaxfm_index); %
    if make_plots
        subplot(2,3,3)
        plot(t,aif_data(ymaxfm_index,:),'k'); hold on;
        xlabel('t [s]'); ylabel('AIF(t)')
        title('max value & lowest fm')
    end
    
    % calculate mean aif
    mean_aif = mean(aif_data(ymaxfm_index,:));
    if make_plots
        subplot(2,3,4)
        plot(t,mean_aif,'r','linewidth',1.5)
        xlabel('t [s]'); ylabel('AIF(t)')
        title('mean aif (max val, min fm)')
    end
catch
    ymaxfm_index = ymax_index; %select only those with max values
end

if make_plots
    subplot(2,3,5)
    plot(t,aif_data_c(ymaxfm_index,:),'k'); hold on;
    xlabel('t [s]'); ylabel('AIF(t)')
    title('curves with highest values')
end

mean_aif_c = mean(aif_data_c(ymaxfm_index,:),1); %mean of the max value curves
if make_plots
    subplot(2,3,6)
    plot(t,mean_aif_c,'r','linewidth',1.5)
    xlabel('t [s]'); ylabel('C(t)')
    title('mean aif')
end
aif_area = trapz(t,mean_aif_c);

if save_aif; save('mean_aif_c.mat', 'mean_aif_c'); end

end