function aif_arrival_time = find_aif_arrival_time(handles)
% THIS IS A REDUDANT FUNCTION. If needed, contact
%              Arthur Simbarashe Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mean_aif_c = handles.mean_aif_c;
tr = handles.tr;
img_size = handles.img_size;

t = 0:tr:(img_size(4)-1)*tr;
dt = t(2)-t(1); 
tl = length(t);

% first we smooth the curve            %edge detection methods are noise-sensitive
fwhm = 4; %filter witdh of kernel1 (smoothing gaussian)
interpolFactor = 2; %linearly interpolate tissue curve and AIF with 2 points between each sample
intFactor = ((tl-1)*interpolFactor + tl); %the interpolation time vector  
elements = 11; %use 7 elements kernels
kernel1 = (0:elements-1)-((elements-1)/2); %-3,-2-1,0,1,2,3
kernel1 = exp(-log(2.0)/(fwhm/2.0)^2.*kernel1.^2); %gaussian filter kernel
%This is just the Gaussian formula without the prefactor before the
%exponential. It comes from that FWHM = 2*sqrt(2*ln2)*sigma where sigma is
%the standard deviation. The mean is zero because our x-array is symmetrical
%about zero. So the expression above is just k = e^(-x^2/2sigma^2). The
%prefactor is not necessary because of the normalisation that happens in
%the next step.
kernel1 = kernel1./sum(kernel1); %normalize
kernel2 = [2.0,0.0,-2.0]; %edge enhancement

aif_edge = interp1(([0.0,0.0,0.0,squeeze(mean_aif_c),0.0,0.0,0.0]),linspace(1,tl+6,intFactor));

aif_edge(aif_edge<0) = 0; %remove negative values 

%perform Canny edge detection
aif_edge = conv(aif_edge, kernel1, 'same');
aif_edge = conv(aif_edge, kernel2, 'same');

%the x-axis corresponding to the interpolated/filtered data
tmpXaxis = [t(1)-3.0*dt,t(1)-2.0*dt,t(1)-dt, t, t(tl)+dt, t(tl)+2.0*dt,t(tl)+3.0*dt];
XaxisInt = linspace(tmpXaxis(1),tmpXaxis(end),intFactor); 

idx = min(find(aif_edge > max(aif_edge)/2,1)); %find where the filtered ncr-curve exeeds 0.33 of its max
%because there we most likely have an edge.
aif_arrival_time = XaxisInt(idx); aif_arrival_time(aif_arrival_time<0) = 0;
end