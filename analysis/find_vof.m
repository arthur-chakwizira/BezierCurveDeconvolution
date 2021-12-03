function [mean_vof_c, vof_area] = find_vof(data, config)

dsc_data = data.dsc_data;
mask = data.mask;
te = config.te; %s
tr =  config.tr; %s
img_size = config.img_size;
baseline_idx = config.baseline_index;
tail_idx = config.tail_index;

vof_config = config.vof_config;
vof_slice = vof_config.vof_slice;
make_plots = vof_config.make_plots;
save_vof = vof_config.save_vof;
vof_time_point = vof_config.vof_time_point;

t = 0:tr:(img_size(4)-1)*tr; %time points at which each of the DSC-MRI images were taken

% start to search for vof:s in manually drawn ROI
figure('position',[200 200 1200 800])
subplot(1,2,1)
colormap(gray) 

% show slice
imagesc(squeeze(dsc_data(:,:,vof_slice,vof_time_point))), axis image off 
%This vof_slice is chosen at the time point corresponding to the minimum on the whole-brain signal curve

h_roi = drawfreehand; %create draggable freehand region
vof_roi = h_roi.createMask; %createMask is a function of the imfreehand object. It creates a mask
%within the image: ones in ROI and zeros outside.
vof_mask = false(img_size(1:3));

% pick out curves from the ROI
vof_slices = (vof_slice-1):(vof_slice+1); %include 3 slices
for i = 1:length(vof_slices) 
    vof_mask(:,:,vof_slices(i)) = vof_roi.*mask(:,:,vof_slices(i));
end

tmp_dsc_data = reshape(dsc_data,[prod(img_size(1:3)) img_size(4)]); 
tmp_mask = repmat(vof_mask(:),[1 img_size(4)]); 
vof_data = reshape(tmp_dsc_data(tmp_mask),[sum(vof_mask(:)) img_size(4)]);

% search for vof based on FM and Ymax

if make_plots
    figure(2)
    plot(t,mean(vof_data),'k'), hold on
    xlabel('t [s]'); ylabel('vof(t)')
    title('mean vof (raw signal)')
end

% convert to concentration
vof_data_c = vof_data.*0; %

% find the indicies with top 5% highest peaks
ymax    = zeros(sum(vof_mask(:)),1);    
fm      = zeros(sum(vof_mask(:)),1);
noise   = zeros(sum(vof_mask(:)),1);
% pointy  = zeros(sum(vof_mask(:)),1);

for i = 1:sum(vof_mask(:))
    
    vof_s0 = mean(mean(vof_data(i,baseline_idx)));

    vof_data_c(i,:) = -(1/te) .* log(vof_data(i,:)./vof_s0) ;
    vof_data(i,:) = vof_data(i,:);
    
    ymax(i) = max(vof_data_c(i,:));
    fm(i) = sum(t.*vof_data_c(i,:))/sum(vof_data_c(i,:));
    
    % noise level
    noise(i) = std(vof_data(i,baseline_idx)) + std(vof_data(i,tail_idx));
    
    % remove to pointy curves
%     pointy(i) = sum(abs(diff(vof_data(i,:))));
end


% high peaks ymax
numvofs = 5;
[~,sortIndex] = sort(ymax(:),'descend');
ymax_index = sortIndex(1:numvofs); % highest values

% remove noisy and pointy
% ymax_index(noise(ymax_index)>40) = [];
% ymax_index(pointy(ymax_index)>2200) = [];

% remove low s0
% ymax_index(vof_s0(ymax_index)<600) = [];


if make_plots
    figure('position',[200 200 1400 400])
    subplot(2,3,1)
    plot(t,vof_data_c(ymax_index,:),'k'), hold on, %ylim([0 120])
    xlabel('t [s]'); ylabel('vof(t)')
    title('curves with highest values')
end


% high fm
numvofs = round(sum(vof_mask(:))/5);
[~,sortIndex] = sort(fm(:),'descend');
fm_index = sortIndex(1:numvofs); 

% remove noisy and pointy
%fm_index(noise(fm_index)>10) = [];
%fm_index(pointy(fm_index)>400) = [];

% remove low s0
% ymax_index(vof_s0(ymax_index)<600) = [];

if make_plots
    subplot(2,3,2)
    plot(t,vof_data_c(fm_index,:),'k'); hold on;
    xlabel('t [s]'); ylabel('vof(t)')
    title('lowest fm values')
end

% find the intersect
try
    ymaxfm_index = intersect(ymax_index,fm_index); %curves with lowest fm and max values
%     numvofs = length(ymaxfm_index); %
    if make_plots
        subplot(2,3,3)
        plot(t,vof_data(ymaxfm_index,:),'k'); hold on;
        xlabel('t [s]'); ylabel('vof(t)')
        title('max value & lowest fm')
    end
    
    % calculate mean vof
    mean_vof = mean(vof_data(ymaxfm_index,:));
    if make_plots
        subplot(2,3,4)
        plot(t,mean_vof,'r','linewidth',1.5)
        xlabel('t [s]'); ylabel('vof(t)')
        title('mean vof (max val, min fm)')
    end
catch
    ymaxfm_index = ymax_index; %select only those with max values
end

if make_plots
    subplot(2,3,5)
    plot(t,vof_data_c(ymaxfm_index,:),'k'); hold on;
    xlabel('t [s]'); ylabel('vof(t)')
    title('curves with highest values')
end

mean_vof_c = mean(vof_data_c(ymaxfm_index,:),1); %mean of the max value curves
if make_plots
    subplot(2,3,6)
    plot(t,mean_vof_c,'r','linewidth',1.5)
    xlabel('t [s]'); ylabel('C(t)')
    title('mean vof')
end
vof_area = trapz(t,mean_vof_c);

if save_vof; save('mean_vof_c.mat', 'mean_vof_c'); end

end