%% main script for simulation of DSC-MRI data
%Author:
%Arthur Chakwizira
%Lund University, Sweden
%arthur.chakwizira@med.lu.se

clc
%% SIMULATE DSC DATA
    cbf_g = [10 20 30 40 50 60 70]; %cbf gray matter
    sim_config.cbv_gray = 0.04;%ml/g
    cbf_w = [5 10 15 20 25 30 35]; %cbf white matter
    sim_config.cbv_white = 0.0200; %ml/g
    sim_config.snr_t = 20; %snr
    sim_config.image_size = [64 64];    %matrix size to simulate
    sim_config.n_slices = 1;            %number of slices
    sim_config.n_time_points = 161;      %number of time points
    sim_config.te =  29E-3;          %TE in seconds
    sim_config.tr = 1.243;          %TR in seconds
    
    %we will try to generate delay and dispersion data automatically
    %[disperse_aif, delay_aif , disp_level, delay_level]
    del_disp = ["false" "false"   "high" "high"; %this line means no delay or dispersion
        
                "false"  "true"   "high"  "high";
                "false" "true"    "high"  "medium";
                "false" "true"    "high"  "low";
    
                "true"  "false"   "high"  "high";
                "true"  "false"   "medium"  "high";
                "true"  "false"  "low"  "high";
    
                "true"  "true"    "high"  "high";
                "true"  "true"    "medium"  "medium";
                "true"  "true"    "low"  "low"];
    
    
for i = 1:length(cbf_g)
    
    sim_config.cbf_gray = cbf_g(i);
    sim_config.cbf_white = cbf_w(i);
    
    for setting = 1:1%size(del_disp, 1)
        settings = del_disp(setting,:);
        sim_config.disperse_aif = eval(settings(1));     %true = include dispersion
        sim_config.delay_aif = eval(settings(2));        %true = include delay
        
        %use only GM for delay or dispersion study
        if (sim_config.disperse_aif||sim_config.delay_aif)
               sim_config.cbv_white = sim_config.cbv_gray; %ml/g
               sim_config.cbf_white =  sim_config.cbf_gray;
        end
        
        %--------------------------------------------------------------
        %select dispersion kernel.
        %choices are 'exponential', 'gamma_dist' and 'lognormal'
        sim_config.dk = 'exponential';
        %select level of dispersion
        %choises are 'high', 'medium' and 'low'
        sim_config.disp_level = settings(3);%'high';
        %-------------------------------------------------------------
        %select level of delay
        %choises are 'high', 'medium' and 'low'
        sim_config.delay_level = settings(4);%'high';
        %-------------------------------------------------------------
        %select residue function
        %choices are 'default' (monoexponential) and 'gamma_dist'
        sim_config.residue_function = 'gamma_dist';
        sim_config.lambda = 1; %shape parameter for gamma residue function
        %-------------------------------------------------------------
        sim_config.make_plots = true;       %display simulated dsc data
        sim_config.save_data = true;        %save simulated dsc data
        %------------------------------------------------------------
        simulate_dsc_data(sim_config)
    end
end

