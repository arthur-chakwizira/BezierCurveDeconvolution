%% main script for analysis of simulated DSC-MRI data using Bezier Curve Deconvolution and oSVD/sSVD
%Specifically for mass deconvolving simulated delay-dispersion data


%Author:
%Arthur Chakwizira
%Department of Medical Radiation Physics
%Lund University, Sweden
%arthur.chakwizira@med.lu.se


close all force
clc
cprintf('*blue', 'ANALYSIS OF DSC-MRI DATA \n')
cprintf('*blue','\n')

%Loop


% %[disp, delay , disp_level, delay_level, disp_correction, delay_correction]

del_disp_corr = ["false" "false"   "high"    "high"  "false"  "false";
    "false" "false"   "high"    "high"  "false"  "true";
    "false" "false"   "high"    "high"  "false"  "true";
    "false" "false"   "high"    "high"  "false"  "true";
    
    "false"  "true"   "high"  "high"   "false"  "false";
    "false"  "true"   "high"  "high"   "false"  "true";
    "false" "true"    "high"  "medium"  "false"  "false";
    "false" "true"    "high"  "medium"  "false"  "true";
    "false" "true"    "high"  "low"     "false"    "false";
    "false" "true"    "high"  "low"     "false"    "true";
    
    "true"  "false"   "high"  "high"   "false"     "false";
    "true"  "false"   "high"  "high"   "false"     "true";
    "true"  "false"   "medium"  "high"   "false"     "false";
    "true"  "false"    "medium"  "high"  "false"     "true";
    "true"   "false"  "low"   "high"    "false"  "false";
    "true"   "false"  "low"   "high"    "false"   "true";
    
    "true"  "true"    "high"  "high"     "false"     "false";
    "true"  "true"    "high"  "high"      "false"     "true";
    "true"  "true"    "medium"  "medium"   "false"   "false";
    "true"  "true"    "medium"  "medium"       "false"   "true";
    "true"  "true"        "low"  "low"    "false"    "false";
    "true"  "true"      "low"  "low"    "false"    "true"];


sn_r = 20;

load('simulated_data/kappa_a.mat')
load('simulated_data/kappa_s.mat')

for cbf_gray = 10:10:70
    %Deconvolve all combinations of delay and dispersion for each cbf, at once
    for setting = 1:size(del_disp_corr, 1)
        settings = del_disp_corr(setting, :);
        
        
        dis = eval(settings(1)); %whether there is dispersion in the input data
        del = eval(settings(2)); %whther there is delay in input data
        
        disp_level =  settings(3);%level of dispersion in data
        del_level = settings(4); %level of delay in data
        
        
        if dis; dis_key = disp_level; else; dis_key = 'no';  end
        if del; del_key = del_level; else; del_key = 'no';  end
        
        
        config.with_dispersion = eval(settings(5));      %true = correct for dispersion
        config.with_delay = eval(settings(6));  % true = correct for delay
        
        if config.with_dispersion && config.with_delay; corr_key = 'both'; end
        if config.with_dispersion && ~config.with_delay; corr_key = 'disp'; end
        if ~config.with_dispersion && config.with_delay; corr_key = 'del'; end
        if ~config.with_dispersion && ~config.with_delay; corr_key = 'none'; end
        
        
        % for run = 1:length(cbf_gray)
        config.target_folder = strcat('SIM-DATA-SNR-', num2str(sn_r) , '/delay_dispersion/cbf_' ,num2str(cbf_gray), '/only_delay_correction/variable_files_deldisp_study_', del_key, '_del_', dis_key,'_disp_', 'corr_', corr_key, '_cbf', num2str(cbf_gray));
        
        path = strcat('SIM-DATA-SNR-', num2str(sn_r) , '/delay_dispersion/cbf_' ,num2str(cbf_gray)); %cbf_' ,num2str(cbf_gray));
        
        file = strcat('St_disp_', num2str(dis), '_exponential_lvl_', disp_level, ...
            '_delay_', num2str(del), '_lvl_', del_level ,'_resfunc_gamma_dist_lambda_1_cbf_', num2str(cbf_gray), '.nii');
        
        path_to_data = fullfile(path, file);
        path2 = strcat('SIM-DATA-SNR-', num2str(sn_r)); file2 = 'aif.nii';
        path_to_aif = fullfile(path2, file2);
        
        %% CONFIGURATION SETTINGS <>
        %Choosing input data
        config.simulate_data = true;       %true = simulate dsc mri data
        config.load_existing_data = true;   %true = load existing dsc mri data
        config.make_mask = true;            %true = generate mask for loaded data
        config.mask_threshold = 0;%100;        %set mask threshold
        config.save_mask = false;            %true = save generated mask
        %Select TE and TR and kappa
        config.te = 29E-3;                  %seconds
        config.tr = 1.243;                  %seconds
        config.kappa = 1;%0.701;%0.705;               %cm^3/g
        %For whole-brain curve, AIF & VOF
        config.create_whole_brain_curve = false;   %true = make whole-brain curve
        config.calculate_aif = false;       %if false, program will request AIF
        config.calculate_arrival_of_aif = false; %AIF arrival time
        config.make_plots = true;          %true = plot whole-brain curve, AIF, VOF
        %Baseline and tail indices
        config.baseline_index = 1:10;
        config.tail_index = 130:161;
        %Select deconvolution algorithm
        config.BzD = true;                  %true = use Bezier curve deconvolution
        config.oSVD = false;                        %true = use oSVD deconvolution
        %Choose whether to correct for delay and/or dispersion
        % config.with_dispersion = false;      %true = correct for dispersion
        % config.with_delay = false;           % true = correct for delay
        %Choose slices to work with, otherwise all available will be calculated
        config.slice_range = [1  10  1  10  1  1];
        %Choose whether to plot residue functions (WARNING: this will take time)
        config.plot_residue_funcs = false;
        %Choose to get notified upon finishing (WARNING: program will play sound)
        config.notify_when_done = true;
        config.save_fitted_params = true; %decide whether to save fitted parameters
        config.save_results = true;     %save results of post-analysis calculations

        
        %% LOAD EXISTING DSC DATA
        if config.load_existing_data
            [dsc_data, mask] = load_dsc_data(path_to_data, config);
            data.dsc_data = dsc_data;
            data.mask = mask;
            config.img_size = size(dsc_data);
        end
        
        %% Target folder
        if ~exist(config.target_folder, 'dir'); mkdir(config.target_folder); end
        
        %% GENERATE WHOLE-BRAIN SIGNAL CURVE
        if config.create_whole_brain_curve
            [mean_dsc_signal, t_min_signal] = whole_brain_curve(data, config);
        end
        
        %% FIND AIF (ROI-BASED METHOD)
        if config.calculate_aif
            aif_config.aif_slice = ceil(config.img_size(3)/2);
            aif_config.make_plots = true;
            aif_config.save_aif = true;
            if ~exist('t_min_signal','var')
                config.make_plots = false;
                [~, t_min_signal] = whole_brain_curve(data, config);
                config.make_plots = true;
            end
            aif_config.aif_time_point = t_min_signal;
            config.aif_config = aif_config;
            [mean_aif_c, aif_area] = find_aif(data, config);
            config.aif_config = [];
        else
            if config.simulate_data
                mean_aif_s = niftiread(path_to_aif);
                aif_s0 = mean(mean_aif_s(config.baseline_index));
                mean_aif_c = -(1/(config.te)) .* log(mean_aif_s./aif_s0);
            else
                mean_aif_c = niftiread(path_to_aif);
            end
        end
        % plot(mean_aif_c)
        data.mean_aif_c = mean_aif_c;

        %% FIND ARRIVAL TIME OF AIF (CANNY EDGE)
        %The Canny edge detector is an operation that uses a multistage
        %algorithm to detect edges in images (places where the signal changes
        %sharply)
        if config.calculate_arrival_of_aif
            aif_arrival_time = find_aif_arrival_time(data,config);
        end
        
        %% BEZIER CURVE DECONVOLUTION
        if config.BzD
            cprintf('*blue','BEZIER CURVE DECONVOLUTION \n')
            cprintf('*blue','\n')
            
            params = bezier_deconvolve(data,config);
            
            fitd_omega     = params(:,:,:,1:5);
            fitd_cbf      	= params(:,:,:,6).*6000; %$$$$%to convert to ml/100g/min
            
            if config.with_dispersion && ~config.with_delay;   fitd_dk  =  params(:,:,:,7:8);  fitd_delay = false;  end
            if config.with_delay && ~config.with_dispersion;   fitd_delay = params(:,:,:,7);  fitd_dk = false;  end
            if config.with_delay && config.with_dispersion
                fitd_delay = params(:,:,:,7);
                fitd_dk  =  params(:,:,:,8:9);
            end
            if ~config.with_delay && ~config.with_dispersion
                fitd_delay  	= false;
                fitd_dk = false;
            end
            
            if config.simulate_data
                kappa = config.kappa;
                %     load('simulated_data/kappa_a.mat')
                %     load('simulated_data/kappa_s.mat')
                %     kappa_a = 2.05;
                % kappa_s = 36;
                fitd_cbf = fitd_cbf.*((kappa_a*kappa)/kappa_s);
            end
            
            niftiwrite(fitd_omega, strcat(config.target_folder, '/fitd_omega_cBzD.nii'))
            niftiwrite(fitd_cbf, strcat(config.target_folder, '/fitd_cbf_cBzd.nii'))
            if config.with_delay; niftiwrite(fitd_delay, strcat(config.target_folder, '/fitd_delay_cBzD.nii')); end
            if config.with_dispersion; niftiwrite(fitd_dk, strcat(config.target_folder, '/fitd_dk_cBzD.nii')); end
        else
            fitd_cbf = false; fitd_delay = false; fitd_dk = false; fitd_omega = false;
        end
        
        %% do oSVD deconvolution
        if config.oSVD
            cprintf('*blue','oSVD DECONVOLUTION \n')
            cprintf('*blue', '\n')
                %svd_outcome = oSVD_deconvolve(dsc_data, mask,mean_aif_c,save_fitted_params, notify_when_done, slice_range);
                if  ~exist('fitd_delay','var') && ~exist('fitd_dk','var')
                    fitd_delay = false; fitd_dk = false;
                end
                [fitd_r_svd, fitd_cbf_svd, fitd_delay_svd] = oSVD_deconvolution(data, config);
                if fitd_r_svd == false & fitd_cbf_svd == false & fitd_delay_svd == false
                    config.oSVD = false;
                end
            
            if config.simulate_data
                kappa = 1; %must be same as in data for oSVD.
                %      load('simulated_data/kappa_a.mat');
                %     load('simulated_data/kappa_s.mat');
                %     kappa_a = 2.05;%0.108;
                %     kappa_s = 36;%27.8;%30;%36;%1.8; %gives 40% signal drop when cbv = 4% and cbf = 60ml/100g/min
                fitd_cbf_svd = fitd_cbf_svd.*((kappa_a)/(kappa_s*kappa)); %for simulated data
            end
         
            
            niftiwrite(fitd_cbf_svd, strcat(config.target_folder, '/fitd_cbf_svd.nii'));
            
        else
            fitd_cbf_svd = false; fitd_r_svd = false; fitd_delay_svd = false;
        end
        
        
        
        
        
        %% plot residue functions
        data.fitd_omega = fitd_omega;
        data.fitd_r_svd = fitd_r_svd;
        
        if config.plot_residue_funcs
            if ~config.oSVD; fitd_r_svd = false; end %if not plotting oSVD, set default value for fitd_r_svd
            if ~exist('fitd_r_svd','var'); fitd_r_svd = false; end
            if ~config.BzD; fitd_omega = false; end
            data.fitd_omega = fitd_omega;
            data.fitd_r_svd = fitd_r_svd;
            plot_residue_functions(data, config)
        end
        
        
        %% calculate CBV and MTT (two ways!), as well as R10/R50 and TTP
        
        if config.simulate_data
            kappa = config.kappa;
            %     kappa_a = 2.05; kappa_s = 36;
            %     load('simulated_data/kappa_a')%kappa_a = 2.05;%0.108;
            %     load('simulated_data/kappa_s')%kappa_s = 30;%27.8;%30;%36; %50;%1.8; %gives 40% signal drop when cbv = 4% and cbf = 60ml/100g/min
            cbv_without_mtt = cbv_no_mtt(data, config).*(kappa_a/(kappa_s*kappa)); %ml/100g kappas for simulated data
            niftiwrite(cbv_without_mtt, strcat(config.target_folder, '/cbv_without_mtt.nii'));
        else
            cbv_without_mtt = cbv_no_mtt(data, config).*cf; %ml/100g
            niftiwrite(cbv_without_mtt, strcat(config.target_folder, '/cbv_without_mtt.nii'));
        end
        
        data.cbv_without_mtt = cbv_without_mtt;
        data.fitd_cbf = fitd_cbf;
        data.fitd_cbf_svd = fitd_cbf_svd;
        
        [mtt_BzD_cbv, mtt_oSVD_cbv] = mtt_with_CBV(data,  config); %seconds
        
        [mtt_BzD_no_cbv, mtt_oSVD_no_cbv] = mtt_no_CBV(data,  config);%seconds
        
        data.mtt_BzD_no_cbv = mtt_BzD_no_cbv;
        data.mtt_oSVD_no_cbv = mtt_oSVD_no_cbv;
        
        [cbv_BzD_mtt, cbv_oSVD_mtt] = cbv_using_mtt(data,  config); %ml/100g
        
        data.fitd_dk = fitd_dk;
        data.fitd_delay = fitd_delay;
        data.fitd_delay_svd = fitd_delay_svd;
        
        [ttp_from_data, ttp_fit_BzD, ttp_fit_oSVD] = find_ttp(data, config); %seconds
        
        [r10_BzD, r50_BzD, r10_svd, r50_svd] = r10_r50(data, config); %in seconds
    end
end


%% calculate OEF and CMRO2
%OEF is oxygen extraction fraction. It reports on the balance between
%oxygen consumption and oxygen delivery. Slight changes in OEF can be a
%sign of physiological perturbation. OEF is defined as
%OEF = oxygen consumption/oxygen delivery = CMRO2/Ca*CBF
%CMRO2 is cerebral metabollic rate of oxygen
save_results = true;
[oef_BzD, cmro2_BzD] = oef_cmro2(fitd_omega,fitd_cbf, img_size, mask, slice_range, save_results);

%% check correlations
% corerlation between delay and ttp (should correlate)
make_plots = true;
slice = 9;
[rho_ttp_delay_BzD, p_ttp_delay_BzD] = correlate(ttp_fit_BzD, fitd_delay, mask, slice, make_plots);
[rho_ttp_delay_oSVD, p_ttp_delay_oSVD] = correlate(ttp_fit_oSVD, fitd_delay_svd, mask, slice, make_plots);

%correlation between mtt and delay?
[rho_mtt_delay_BzD, p_mtt_delay_BzD] = correlate(mtt_BzD_no_cbv, fitd_delay, mask, slice, make_plots);
[rho_mtt_delay_oSVD, p_mtt_delay_oSVD] = correlate(mtt_oSVD_no_cbv, fitd_delay_svd, mask, slice, make_plots);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%