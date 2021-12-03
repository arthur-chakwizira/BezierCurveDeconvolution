function simulate_dsc_data(sim_config)
%% Generating simulated DSC-MRI data
%Author:
%Arthur Chakwizira
%Department of Medical Radiation Physics
%Lund University, Sweden
%arthur.chakwizira@med.lu.se


cprintf('*blue','SIMULATING DSC-MRI DATA. \n')
cprintf('*blue','\n')

%retrieve needed data
%--------------------------------------------------------------------
te = sim_config.te; %s
tr = sim_config.tr; %s
image_size = sim_config.image_size;
n_slices = sim_config.n_slices;
n_time_points = sim_config.n_time_points;
delay_aif = sim_config.delay_aif;
delay_level = sim_config.delay_level;
disperse_aif = sim_config.disperse_aif;

dk = sim_config.dk;
disp_level = sim_config.disp_level;
residue_function = sim_config.residue_function;
lambda = sim_config.lambda;
make_plots = sim_config.make_plots;
save_data = sim_config.save_data;

snr_t = sim_config.snr_t;

cbv_gray = sim_config.cbv_gray;%ml/g
cbv_white = sim_config.cbv_white; %ml/g
cbf_gray = sim_config.cbf_gray/6000; %ml/g/s
cbf_white = sim_config.cbf_white/6000; %ml/g/s
mtt_gray = cbv_gray/cbf_gray;
mtt_white = cbv_white/cbf_white;

kappa = 1; %cm^3/g
%-------------------------------------------------------------------

im_size = [image_size(1), image_size(2),  n_slices, n_time_points];
treal = 0:tr:(im_size(4) - 1)*tr;


t = treal;%0:10E-3:treal(end);
dt = t(2) - t(1);
%aif
t0 = 20; %s
b = 1.5;% s
r = 3.0;
C_0 = 1;
C_a = zeros(1, length(t));
for i = 1:length(t)
    if t(i) <= t0
        C_a(i) = 0;
    else
        C_a(i) = C_0*((t(i) - t0)^r)*exp(-(t(i)-t0)/b);
    end
end

%delay and/or disperse aif
if delay_aif
    switch delay_level
        case 'low'
            delay = 1;
        case 'medium'
            delay = 3;
        case 'high'
            delay = 6;
    end
    tail_idx = im_size(4)-5:im_size(4);
    aif_tail_c = mean(C_a(tail_idx));
    
    C_a_delay = interp1(t, C_a, t-delay, 'pchip', 'extrap');
    if (delay < 0); id = ceil(-delay/dt);
        C_a_delay(end-id:end) = aif_tail_c;
    end
    if (delay > 0); id = ceil(delay/dt);
        C_a_delay(1:id) = 0;
    end
else
    C_a_delay = C_a;
end

%% Convert aif to signal
% kappa_a = -log(0.4)/(te*max(C_a(:)))%
if dt == 10e-3; kappa_a = 6.9644; end %with higher sampling rate (10 ms)
if dt == tr; kappa_a = 7.0281; end %with lower sampling rate (tr)

S_a_0 = 100;
S_a = S_a_0.*exp(-kappa_a.*C_a*te);
S_a = interp1(t, S_a, treal, 'pchip'); %resample to tr
% plot(treal, S_a)

% and add Rice noise
reps = 1;
snr_a = snr_t/sqrt(reps);
sd_noise_a = S_a_0/snr_a;

realchannel_a = zeros(reps, im_size(4));
imaginarychannel_a = zeros(reps, im_size(4));
% S_a_noisy = zeros(reps, im_size(4));

cprintf('*blue','\n')
cprintf('*blue','AIF: \n')
cprintf('*blue','\n')
parfor_progress(reps);
for rep = 1:reps
    realchannel_a(rep, :) = normrnd(0, sd_noise_a, 1, im_size(4)) + S_a;
    imaginarychannel_a(rep, :) = normrnd(0, sd_noise_a, 1, im_size(4));
    parfor_progress;
end
parfor_progress(0);
realnoise = mean(realchannel_a, 1);
imaginarynoise = mean(imaginarychannel_a, 1);
S_a_noisy = sqrt( realnoise.^2 + imaginarynoise.^2  );

% S_a_noisy = mean(S_a_noisy, 1);

%plot aif
if make_plots
    figure()
    hold on
    subplot(3,3,1)
    plot(t,C_a)
    ylabel('[a.u]')
    xlabel('time [s]')
    title('AIF concentration')
    
    subplot(3,3,2)
    plot(treal, S_a)
    ylabel('[a.u]')
    xlabel('time [s]')
    title('AIF signal')
    
    subplot(3,3,3)
    plot(treal, S_a_noisy)
    ylabel('[a.u]')
    xlabel('time [s]')
    title('AIF signal + noise')
end

if save_data
    if ~exist('simulated_data', 'dir')
        mkdir simulated_data
    end
    save('simulated_data/kappa_a.mat', 'kappa_a')
    if ~exist('simulated_data/aif.nii', 'file')
        niftiwrite(S_a_noisy, 'simulated_data/aif.nii')
    end
end


%%
if disperse_aif
    switch dk
        case 'exponential'
            switch disp_level
                case 'low'
                    betaa = 1/1.5; %/s. 1/beta is the effective mtt from the site of aif measurement to the input to the particular voxel
                case 'medium'
                    betaa = 1/3.0;
                case 'high'
                    betaa = 1/4.5;
            end
            disp_kernel = betaa.*exp(-betaa.*t);
            C_a_delay_disp = dt.*filter(C_a_delay, 1, disp_kernel); %we will use conv instead of filter.
            plot(t, C_a_delay, t,C_a_delay_disp, 'g-')
            
        case 'lognormal'
            switch disp_level
                case 'low'
                    sigma = 1; %s. This is from  Mendirahta modeling and correction of bolus dispersion effects in dsc mri
                    mu = -1;
                case 'medium'
                    sigma = 0.75;  mu = -0.15;
                case 'high'
                    sigma = 0.78;  mu = 0.59;
            end
            disp_kernel = (1./(sigma.*t.*sqrt(2*pi))).*exp(-( (log(t) - mu).^2 )./(2*sigma^2));
            disp_kernel(1) = interp1(t(2:end),disp_kernel(2:end),0,'linear','extrap');
            disp_kernel(disp_kernel<0) = 0;
            disp_kernel = disp_kernel./trapz(t,disp_kernel); % normalize
            disp_kernel = [zeros(1,161), disp_kernel];
            C_a_delay_disp = dt.*filter(C_a_delay, 1, disp_kernel); %we will use conv instead of filter.
            whos C_a_delay_disp
        case 'gamma_dist'
            switch disp_level
                case 'low'
                    p = 1; s = 2; %from same paper as lognormal parameters
                case 'medium'
                    p = 3; s = 1;
                case 'high'
                    p = 5; s = 0.5;
            end
            disp_kernel = ((s^(1+s*p))/gamma(1+s*p)).*t.^(s*p).*exp(-s.*t);
            C_a_delay_disp = dt.*filter(C_a_delay, 1, disp_kernel); %we will use conv instead of filter.
    end
    
    %     figure; hold on; plot(C_a_delay_disp, 'r-'); plot( C_a, 'b-')
    
else
    C_a_delay_disp = C_a_delay; %don't disperse aif if told  not to
end


%tissue concentration curves and residue function

switch residue_function
    case 'default'
        R_gray = exp(-t./mtt_gray);
        R_white = exp(-t./mtt_white);
    case 'monoexponential'
        t_fine = 0:0.001:t(end);
        %         t_fine = t;
        R_gray = exp(-t_fine./mtt_gray);
        R_white = exp(-t_fine./mtt_white);
        R_gray = interp1(t_fine, R_gray,t, 'pchip');
        R_white = interp1(t_fine, R_white,t, 'pchip');
    case 'gamma_dist'
        t_fine = 0:0.01:treal(end);
        %           t_fine = t;
        beeta_gray = mtt_gray/lambda;
        beeta_white = mtt_white/lambda;
        h_gray = @(time) (1/( (beeta_gray^lambda)*gamma(lambda) )).*(time.^(lambda-1)).*exp(-time./beeta_gray);
        h_white = @(time)  (1/( (beeta_white^lambda)*gamma(lambda) )).*(time.^(lambda-1)).*exp(-time./beeta_white);
        R_gray = zeros(1, length(t_fine));
        R_white = zeros(1, length(t_fine));
        k = 10;
        for i = 1:length(t_fine)
            if lambda > 5; k = 6; end
            R_gray(i) =  integral(h_gray, t_fine(i), k*t_fine(end));%trapz(t_fine(i:end), h_gray(i:end));  trapz(th(i:end), h_gray(i:end));
            R_white(i) =  integral(h_white, t_fine(i), k*t_fine(end));%trapz(t_fine(i:end), h_white(i:end));   trapz(th(i:end), h_white(i:end));
            %             dbstop if infnan
        end
        if make_plots
            subplot(3,3,4)
            hold on
            plot(t_fine, R_gray,'r-', 'LineWidth', 2)
            
            plot(t_fine, R_white, 'k-', 'LineWidth', 2)
            xlim([0 50])
            %         ylabel('[a.u]')
            xlabel('time [s]')
            title('RESIDUE FUNCTIONS')
            legend('Gray matter', 'White matter')
            hold off
        end
        R_gray = interp1(t_fine, R_gray,t, 'pchip');
        R_white = interp1(t_fine, R_white,t, 'pchip');
        %This model is from Mouridsen et al. Bayesian estimation of cerebral perfusion using a physiological
        %model of microvasculature
        
    case 'linear'
        t_fine = 0:0.1:t(end);
        %          t_fine = t;
        R_gray = zeros(1, length(t_fine));
        R_white = zeros(1, length(t_fine));
        for i = 1:length(t_fine)
            if t_fine(i) <= 2*mtt_gray
                R_gray(i) = 1 - (t_fine(i)/(2*mtt_gray));
            else
                R_gray(i) = 0;
            end
            if t_fine(i) <= 2*mtt_white
                R_white(i) = 1 - (t_fine(i)/(2*mtt_white));
            else
                R_white(i) = 0;
            end
        end
        if make_plots
            %         figure(2)
            subplot(3,3,4)
            hold on
            plot(t_fine, R_gray,'k-')
            plot(t_fine, R_white, 'r-')
            xlim([0 30])
            title('RESIDUE FUNCTIONS')
            legend('Gray matter', 'White matter')
        end
        hold off
        R_gray = interp1(t_fine, R_gray,t, 'pchip');
        R_white = interp1(t_fine, R_white,t, 'pchip');
end


%initialising tissue concentration arrays
C_t_gray = zeros(im_size(1)/2, im_size(2)/2, length(t));
C_t_white = zeros(im_size(1)/2, im_size(2)/2, length(t));
C_t_gray_del_disp = zeros(im_size(1)/2, im_size(2)/2, length(t));
C_t_white_del_disp = zeros(im_size(1)/2, im_size(2)/2, length(t));


for x = 1:im_size(1)/2
    for y = 1:im_size(2)/2
%                 C_t_gray(x,y,:) = kappa*cbf_gray*dt*filter(C_a, 1, R_gray);
%                 C_t_white(x,y,:) = kappa*cbf_white*dt*filter(C_a, 1, R_white);
        
        C_t_gray_del_disp(x,y,:) = kappa*cbf_gray*dt.*filter(C_a_delay_disp,1,  R_gray);
        C_t_white_del_disp(x,y,:) = kappa*cbf_white*dt.*filter(C_a_delay_disp, 1, R_white);
        
        C_t_gray(x,y,:) = C_t_gray_del_disp(x,y,:);
        C_t_white(x,y,:) = C_t_white_del_disp(x,y,:);
        
    end
end


%% create tissue signal matrix
S_0 = 100;
% kappa_s = -log(0.6)/(te*max(C_t_gray(:))) %gives 40% signal drop when cbv = 4% and cbf = 60ml/100g/min

% higher sampling rate (10 ms)
if dt == 10e-3
    switch lambda
        case 100; kappa_s = 107.2732;
        case 5; kappa_s = 118.9152;
        case 1; kappa_s = 151.3207;
    end
% kappa_s = 107.2732; %lambda = 100
% kappa_s = 118.9152; % lambda = 5
% kappa_s = 151.3207; % lambda = 1
end


%lower sampling rate (tr)
if dt == tr
    switch lambda
        case 100; kappa_s = 96.5102;
        case 5; kappa_s = 106.8595;
        case 1; kappa_s = 130.8347;
    end
% kappa_s = 110.7134; %linear
% kappa_s = 96.5102; %lambda = 100
% kappa_s = 106.8595; %lambda = 5
% kappa_s = 130.8347; %lambda = 1
end


% S_t = zeros(im_size);
S_t(:,:,1,:) = S_0.*[exp(-kappa_s.*C_t_gray*te)  , exp(-kappa_s.*C_t_gray_del_disp*te);
    exp(-kappa_s.*C_t_white*te)   ,   exp(-kappa_s.*C_t_white_del_disp*te)];

% S_t(:,:,1,:) = S_0.*[exp(-kappa_s.*C_t_gray*te); exp(-kappa_s.*C_t_gray_del_disp*te)];
S_t_new = zeros(im_size);

%resample to tr
for x = 1:im_size(1)
    for y = 1:im_size(2)
        S_t_new(x,y,1,:) = interp1(t, squeeze(S_t(x,y,1,:)), treal, 'pchip');
    end
end
S_t = S_t_new;

%add Rice noise
sd_noise_t = S_0/snr_t;


S_t_noisy = S_t*0;


cprintf('*blue','C(t): \n')
cprintf('*blue','\n')
parfor_progress(im_size(1));
for x = 1:im_size(1)
    for y = 1:im_size(2)
        realchannel_t = normrnd(0, sd_noise_t, 1, im_size(4)) + squeeze(S_t(x,y,1,:))';
        imaginarychannel_t = normrnd(0, sd_noise_t, 1, im_size(4));
        S_t_noisy(x,y,1,:) = sqrt(realchannel_t.^2 + imaginarychannel_t.^2);
    end
    parfor_progress;
end
parfor_progress(0);



if make_plots
    subplot(3,3,5)
    plot(treal, squeeze(S_t(1,1,1,:)))
    ylabel('[a.u]')
    xlabel('time [s]')
    title('C(t) no noise')
    
    subplot(3,3,6)
    plot(treal,squeeze((S_t_noisy(1,1,1,:))))
    ylabel('[a.u]')
    xlabel('time [s]')
    title('C(t) + noise')
    
    subplot(3,3,7)
    imagesc(squeeze(S_t_noisy(:,:,3)))
    axis image off
    colormap(gray)
    set(gca, 'CLim', [0 120])
    title('Image at 2.5 s')
    
    subplot(3,3,8)
    imagesc(squeeze(S_t_noisy(:,:,22)))
    axis image off
    colormap(gray)
    set(gca, 'CLim', [0 120])
    title('Image at 27 s')
    
    subplot(3,3,9)
    imagesc(squeeze(S_t_noisy(:,:,103)))
    axis image off
    colormap(gray)
    set(gca, 'CLim', [0 120])
    title('Image at 126 s')
    
    
    sgtitle('SIMULATED DSC-MRI DATA')
end

if save_data
    save('simulated_data/kappa_s.mat', 'kappa_s')
    if strcmp(residue_function,'gamma_dist')
        niftiwrite(S_t_noisy, strcat('simulated_data/St_disp_', num2str(disperse_aif), '_', dk, '_lvl_', disp_level, '_delay_', ...
            num2str(delay_aif), '_lvl_', delay_level, '_resfunc_' ,residue_function, '_lambda_', num2str(lambda), ...
            '_cbf_', num2str(cbf_gray*6000), '.nii'))
    elseif strcmp(residue_function,'linear')
        niftiwrite(S_t_noisy, strcat('simulated_data/St_disp_', num2str(disperse_aif), '_', dk, '_lvl_', disp_level, '_delay_', ...
            num2str(delay_aif), '_lvl_', delay_level, '_resfunc_' ,residue_function, ...
            '_cbf_', num2str(cbf_gray*6000), '.nii'))
    else
        niftiwrite(S_t_noisy, strcat('simulated_data/St_', 'disp_', num2str(disperse_aif), '_' ,dk, '_lvl_', disp_level, '_delay_' ,...
            num2str(delay_aif), '_lvl_', delay_level, '_cbf_gw_', num2str(cbf_gray*6000), '_', num2str(cbf_white*6000), '.nii'))
    end
end
cprintf('magenta','Done \n')
cprintf('blue','_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ \n')

end
