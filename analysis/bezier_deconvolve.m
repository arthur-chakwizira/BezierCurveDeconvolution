function params = bezier_deconvolve(data, config)

%This function handles Bezier curve deconvolution

dsc_data = data.dsc_data;
mask = data.mask;
mean_aif_c = data.mean_aif_c;

with_dispersion = config.with_dispersion;
with_delay = config.with_delay;
save_fitted_params = config.save_fitted_params;
notify_when_done = config.notify_when_done;
slice_range = config.slice_range;

te = config.te;
tr = config.tr;
img_size = config.img_size;
target_folder = config.target_folder;

tmp_baseline_idx = config.baseline_index;
tail_idx = config.tail_index;


%if no slice range has been specified, calculate all 20 slices
if  isempty(slice_range)
    slice_range = [1 img_size(3)];
end

t = 0:tr:(img_size(4)-1)*tr;

para.kappa = config.kappa;
para.c_aif = mean_aif_c;
para.aif_tail_c = mean(mean_aif_c(tail_idx));
para.te = te;

gdk_0  = [  2   2 ]; % s p  large p and small s means more dispersion
dk_0  = log(gdk_0); % initialise dk
% dk_0  = (gdk_0); % initialise dk 
delay_0  = 5;


%--------------------------------------------------------------------------
if with_dispersion && ~with_delay  %only dispersion correction
    %for signal modelling
    %     omega_0       = [     15   0.8   8  0.0   25   ];
    %     omega_0 = [10  0.2  10  0.2  100];
    %     omega_0       = [     8   0.5  8  0.5  2  0.0   15   ];
    omega_0       = [     8   0.5   2  0.0   15   ];
    %     omega_0 = [2  0.2  2  0.2  100];
    %     omega_0       = [  3  0.5   2  0.0  100   ];
    f_0 = 1;
    num_pars = 8;
    delay_0 = [];
end
%--------------------------------------------------------------------------
if with_delay && ~with_dispersion   %only delay correction
    
    % omega_0       = [     5   0.5  2  0.2  2  0.0   30   ];
    % omega_0       = [     2   0.9   0.1  0.0   15   ];
    omega_0       = [     8   0.5   2  0.0   15  ];
    %for signal modelling
    %     omega_0       = [     20   0.8   8  0.0   100   ];
    %     omega_0 = [2  0.2  2  0.2  100];
    %     omega_0       = [    3  0.5   2  0.0  100   ];
    f_0 = 0.01;
    num_pars = 7;
    dk_0 = [];
end
%--------------------------------------------------------------------------
if  with_dispersion && with_delay  %delay corr. and dispersion corr.
    %for signal
    %         omega_0 = [10  0.2  10  0.2  100];
    %     omega_0       = [     8   0.5  8  0.5  2  0.0   15   ];
    omega_0       = [     8   0.5   2  0.0   15  ];
    %     omega_0 = [2  0.2  2  0.2  100];
    %     omega_0       = [     3  0.5   2  0.0  100   ];
    f_0 = 1;%10;%0.1;
    num_pars = 9;
end
%--------------------------------------------------------------------------
if ~with_dispersion && ~with_delay  %no delay corr., no dispersion corr.
    %     omega_0       = [       5   0.8   2  0.2   20   ];
    %     omega_0 = [4  0.5  8  0.5  20];
    %     omega_0       = [     8   0.5  8  0.5  2  0.0   15   ];
    %for signal modelling
    %     omega_0       = [     25   0.8   18  0.8   55   ];
    %     omega_0 = [10  0.2  10  0.2  100];
    %     omega_0 = [6  0.5  8  0.2  10];
    % omega_0       = [     15  0.99   1  0.0  15   ];
    omega_0       = [     8   0.5   2  0.0   15 ];
    %  omega_0       = [  0.01   2   0.005   2  0.0   15  0.0  ];
    %     omega_0 = [2  0.2  2  0.2  100 0];
    f_0 = 0.01;
    dk_0 = [];
    delay_0 = []; %We dont want to correct for delay if not correcting for dispersion
    num_pars = 6;
end
%--------------------------------------------------------------------------

x0 = [ omega_0  f_0  delay_0  dk_0  ];
x0(x0==0) = 1E-10;

p0 = log(x0); %estimating parameters on logarithmic scale

xdata = t';

%Default tolerance is 1e-4
% options = optimset('LargeScale', 'off','display','on','TolFun',1E-9,'TolX',1E-9 , 'MaxFunEvals', 20000, 'MaxIter', 20000);
% options = optimoptions(@fminunc, 'OptimalityTolerance', 1E-9);
options = optimset('display','off');

img_size_23 = img_size(2:3);

params = zeros([img_size(1:3) num_pars]);

xrange = slice_range(1):slice_range(2);
yrange = slice_range(3):slice_range(4);
zrange = slice_range(5):slice_range(6);

parfor_progress(img_size(1)); %works together with parfor. Initialises the progress monitor for a
%set of upcoming calculations

for x = xrange%1:img_size(1)
    
    tmp_p = zeros([img_size_23  num_pars]);
    x0 = [ omega_0 f_0  delay_0  dk_0  ];
    x0(x0==0) = 1E-10;
    p0 = log(x0);
    
    for y = yrange%1:img_size_2
        for z = zrange%img_size_3
            if mask(x,y,z)
                
                tmp_tissue_signal   = squeeze(dsc_data(x,y,z,:));
                ydata = tmp_tissue_signal;
                
                tmp_noise_level     = std(ydata(tmp_baseline_idx)); %noise level for this voxel is std of baseline concentration
                
                
                % normalize
                Pnorm = max(p0(:));
                
                % variational bayesian analysis (VBA) fit
                anony_fun = @(p)bezier_fit_function(p,xdata,ydata,para,Pnorm, tmp_noise_level,x0, with_delay, with_dispersion);
                [tmp_tmp_p, ~] = fminunc(anony_fun,p0./Pnorm, options); %normalise initial conditions
                
                tmp_p(y,z,:)     = exp(Pnorm.*tmp_tmp_p); %"undo" the normalisation for further calculations
%                 tmp_p(y,z,:)     = (Pnorm.*tmp_tmp_p); %"undo" the normalisation for further calculations
            end
        end
    end
    
    params(x,:,:,:) = tmp_p;
    
    parfor_progress;
end
parfor_progress(0);

if save_fitted_params
    niftiwrite(params, strcat(target_folder, '/params_cBzD.nii'))
end

if notify_when_done
    beep; pause(0.5); beep; pause(0.5); beep;
end
cprintf('magenta','Bezier curve deconvolution complete. \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')
end