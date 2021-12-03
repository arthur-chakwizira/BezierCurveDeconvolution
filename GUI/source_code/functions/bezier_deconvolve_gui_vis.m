function params = bezier_deconvolve_gui_vis(handles)
%        This function handles Bezier-curve deconvolution. It displays
%        results during computation. Its input is the handles structure and
%        its output is the array params which contains the following fitted
%        parameters: control points for residue functions, cbf, delay and
%        or dispersion kernel.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%retrieve needed data
dsc_data = handles.dsc_data;
mask = handles.mask;
mean_aif_c = handles.mean_aif_c;

with_dispersion = handles.with_dispersion;
with_delay = handles.with_delay;
slice_range = handles.slice_range;

tr = handles.tr;
img_size = handles.img_size;

tmp_baseline_idx = handles.baseline_index;
tail_idx = handles.tail_index;
%__________________________________________________________________________
t = 0:tr:(img_size(4)-1)*tr; 

para.kappa = handles.kappa; 
para.c_aif = mean_aif_c; 
para.aif_tail_c = mean(mean_aif_c(tail_idx));

%Initialise parameters
gdk_0  = [  2   2 ]; % s p  large p and small s means more dispersion
dk_0  = log(gdk_0);
delay_0  = 5;

%--------------------------------------------------------------------------
if with_dispersion && ~with_delay  %only dispersion correction
    omega_0       = [     8   0.5   2  0.0   15   ]; %control points
    f_0 = 1; %cbf
    num_pars = 8; %number of parameters
    delay_0 = []; %delay
end
%--------------------------------------------------------------------------
if with_delay && ~with_dispersion   %only delay correction
    omega_0       = [     8   0.5   2  0.0   15   ];    
    f_0 = 0.01;
    num_pars = 7;
    dk_0 = [];
end
%--------------------------------------------------------------------------
if  with_dispersion && with_delay  %delay corr. and dispersion corr.
    omega_0       = [     8   0.5   2  0.0   15   ];
    f_0 = 1;
    num_pars = 9;
end
%--------------------------------------------------------------------------

if ~with_dispersion && ~with_delay  %no delay corr., no dispersion corr.
    omega_0       = [     8   0.5   2  0.0   15   ];    
    f_0 = 0.01;
    dk_0 = [];
    delay_0 = [];
    num_pars = 6;
end    
%--------------------------------------------------------------------------
x0 = [ omega_0  f_0  delay_0  dk_0  ];  
x0(x0==0) = 1E-10;
p0 = log(x0);
xdata = t';  
%Default tolerance is 1e-4
% options = optimset('LargeScale', 'off','display','off','TolFun',1E-9,'TolX',1E-9 , 'MaxFunEvals', 20000, 'MaxIter', 20000);
% options = optimoptions(@fminunc, 'OptimalityTolerance', 1E-9);
options = optimset('display','off');

img_size_12 = img_size(1:2); 

params = zeros([img_size(1:3) num_pars]);

xrange = slice_range(1):slice_range(2);
yrange = slice_range(3):slice_range(4);
zrange = slice_range(5):slice_range(6);

% Ntot = zrange(end)-zrange(1); 
% if Ntot == 0; Ntot = 1;end
% parfor_progress_gui(handles, 'bzd', Ntot); %initialise progress bar

for z = zrange
    tmp_p = zeros([img_size_12  num_pars]); 
    x0 = [ omega_0  f_0  delay_0  dk_0  ];
    x0(x0==0) = 1E-10;
    p0 = log(x0);
    for x = xrange
        for y = yrange
            if mask(x,y,z) 
                tmp_tissue_signal   = squeeze(dsc_data(x,y,z,:)); 
                ydata = tmp_tissue_signal;
                tmp_noise_level     = std(ydata(tmp_baseline_idx));     
                Pnorm = max(p0(:)); %normalise input for better performance
                try  
                anony_fun = @(p)bezier_fit_function_gui(p,xdata,ydata,para,Pnorm, tmp_noise_level,x0, handles);  %goal function 
                [tmp_tmp_p, ~] = fminunc(anony_fun,p0./Pnorm, options); %solve with fminunc
                tmp_p(x,y,:)     = exp(Pnorm.*tmp_tmp_p);               
                catch
                   tmp_p(x,y,:) = zeros(1,num_pars);
                end
            end
        end
        if handles.show_cbf; display_cbf(tmp_p(x,:,:),x,z,'bzd', handles); end
        if handles.show_mtt; display_mtt(tmp_p(x,:,:),x,z,'bzd', handles); end
        if handles.show_delay; display_delay(tmp_p(x,:,:),x,z,'bzd', handles); end
        if handles.plot_residue_funcs; display_rt(tmp_p(x,:,:),x,z,'bzd', handles); end
    end    
    params(:,:,z,:) = tmp_p;
    
%     parfor_progress_gui(handles, 'bzd');
end
% parfor_progress_gui(handles, 'bzd',0);


wrap_text = 'Bezier-curve deconvolution complete';
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
end