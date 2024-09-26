function [oef_BzD, cmro2_BzD, oef_svd, cmro2_svd] = oef_cmro2_gui(handles)
%       This function accepts the handles structure and returns OEF and
%       CMRO2, for BzD and SVD (using SVD here is not recommended)
%       Author:
%              Arthur Chakwizira
%              arthur.chakwizira@med.lu.se
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


hObject = handles.hObject;
handles = guidata(hObject);

%get needed data__________________________________________________________
fitd_omega = handles.fitd_omega;
fitd_cbf = handles.fitd_cbf;
fitd_cbf_svd = handles.fitd_cbf_svd;
img_size = handles.img_size;
slice_range = handles.slice_range;
target_folder = handles.target_folder;
fitd_r_svd = handles.fitd_r_svd;
tr = handles.tr;
t = 0:tr:(img_size(4)-1)*tr;
%From Mouridsen 2014: Reliable estimation of capillary transit time
%distributions using DSC-MRI
%dC/dx = -k*tau(a_h*P_50(C/(B-C))^(1/h) - a_h*C_t);
%C is the oxygen concentration and the oxygen extraction fraction for a
%given transit time is given by Q = 1-C(1)/C(0)
%Maximum oxygen extraction is then given by
%OEFmax = integral(0, infinity) h(tau)*Q(tau)dtau
%where h(tau) is the distribution of transit tiems (h(t = )-dR/dt)

%ask user to either load OEF calibration or perform calibration
% load_oef_calib = 0;
opts.Interpreter = 'tex';
opts.Default = 'Load calibration';
response = questdlg('\fontsize{10}OEF calculation requires calibration of the rate constant to give desired OEF in healthy white matter. Select action.', ...
    'OEF calibration','Load calibration', 'Perform new calibration', opts);
switch response
    case 'Load calibration'
        load_oef_calib = true;
    case 'Perform new calibration'
        load_oef_calib = false;
        wrap_text = 'Calibrating OEF calculation...';
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
        blink_this(handles.edit7, 'b') %make the text blink
    otherwise
        report('Analysis terminated by user.', handles, 'r')
        oef_BzD   = zeros(img_size(1:3));
        cmro2_BzD   = zeros(img_size(1:3));
        oef_svd   = zeros(img_size(1:3));
        cmro2_svd   = zeros(img_size(1:3));
        return
end

if load_oef_calib
    [file, path] = uigetfile('*.*','Select OEF calibration file', handles.file_folder);
    if file == 0 %if user cancels file selection dialog
        wrap_text = 'No OEF calibration file selected. Program will now request a new calibration.'; %report this on display panel
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
        load_oef_calib = false;
    else
        path_to_data = fullfile(path, file); %otherwise generate full file path
        try load(path_to_data, "k", "Pt", "P50", "h", "OEF_WM"); catch; end %Pt and P50 are also saved to calibration file
        if ~exist('k', 'var') || ~exist('Pt', 'var') || ~exist('P50', 'var') || ~exist('h', 'var')
            report('Selected file is missing k or Pt or P50 or h. New calibration required.', handles, 'r')
            load_oef_calib = false;
            k = NaN;
            Pt = NaN;
            P50 = NaN;
            h = NaN;
        else
            fileID = fopen(strcat(target_folder, '/OEF_calibration_log.txt'), 'w');
            fprintf(fileID, "OEF calibration settings: OEF (healthy white matter) = " + num2str(OEF_WM) + " percent , Pt = " + num2str(Pt) +  " mmHg and P50 = " + num2str(P50) + " mmHg, h = " + num2str(h) + "\n");
            fprintf(fileID, "The calibration above gave k = " + num2str(k) + " /s" + "\n");
            fclose(fileID);
        end
    end
    tmp = get(handles.edit14, 'String');
    tmp{end+1} = '- - - - OEF calibration file - - - -';
    if exist('path_to_data', 'var'); tmp{end+1} = path_to_data; end
    set(handles.edit14, 'String',tmp)
end

if ~load_oef_calib
    
    if ~handles.dsc_data_loaded
        wrap_text = 'OEF calibration requires DSC data. Load data to proceed.'; %check that data exists
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r') %inform user
        return
    end
    
    %ask user for Pt here (Pt is tissue oxygen tension in mmHg) and P50
    %(P50 is oxygen pressure required for half saturation) and h (Hill coefficient)
    wrap_text = 'Waiting for user input...';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    opts.Interpreter = 'tex';
    prompt = {'\color{blue} \fontsize{10} Desired OEF in healthy white matter [%] :', '\color{blue} \fontsize{10} Tissue oxygen tension, Pt [mmHg] :'...
        , '\color{blue} \fontsize{10} Oxygen pressure required for half saturation, P50 [mmHg] :', '\color{blue} \fontsize{10} Hill coefficient, h :'};
    dlgtitle = 'DECONVOLVER: OEF calibration settings';
    dims = [1 90];
    defaults = {'30', '25', '26', '2.8'};
    response = inputdlg(prompt, dlgtitle, dims, defaults, opts);
    if isempty(response) %if use clicks cancel, report and terminate
        wrap_text = 'Target WM OEF, Pt, P50 and h not specified. Proceeding with default values...';
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
        OEF_WM = 30; Pt = 25; P50 = 26; h = 2.8;
    else
        OEF_WM = str2double(response{1});
        Pt = str2double(response{2});
        P50 = str2double(response{3});
        h = str2double(response{4});
        if isfinite(OEF_WM)&&isfinite(Pt)&&isfinite(P50)&&isfinite(h)
            wrap_text = 'User-supplied OEF calibration settings saved. Draw ROI in healthy WM to proceed';
            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
        else
            wrap_text = 'Target WM OEF, Pt, P50 and h not specified.. Proceeding with default values...';
            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
            OEF_WM = 30; Pt = 25; P50 = 26; h = 2.8;
        end
    end
    
    fileID = fopen(strcat(target_folder, '/OEF_calibration_log.txt'), 'w');
    fprintf(fileID, "OEF calibration settings: \n");
    fprintf(fileID,  "OEF (healthy white matter) = " + num2str(OEF_WM) + " percent , Pt = " + num2str(Pt) +  " mmHg, P50 = " + num2str(P50) + " mmHg and h = " + num2str(h) +  "\n");
    fclose(fileID);
    
    
    figure(handles.figure1) %bring main window to top
    
    childlist = get(handles.figure1, 'Children');
    for c = 1:numel(childlist)
        try set(childlist(c), 'Enable', 'on'); catch; end
    end
    
    %% Draw ROIs
    set(handles.pushbutton30, 'String', 'Continue');
    set(handles.pushbutton30, 'Visible', 'on');
    set(handles.pushbutton30, 'Callback', {@Proceed})
    set(handles.pushbutton26, 'Visible', 'on');
    set(handles.pushbutton24, 'Visible', 'on');
    set(handles.pushbutton26, 'Callback', {@addROI});
    set(handles.pushbutton24, 'Callback', {@clearROIs});
    set(handles.slider1, 'Callback', {@change_slice});
    set(handles.slider8, 'Callback', {@change_timepoint});
    prepare_dsc_data(hObject)
    prepare_sliders(hObject)
    initialise_roi_placement(hObject)
    waitfor(handles.pushbutton30, 'Visible')
    handles = guidata(hObject);
    active_rois = handles.rois(ishandle(handles.rois));
    active_roi_slices = handles.roi_slices(ishandle(handles.rois));
    oef_mask = false(handles.img_size(1), handles.img_size(2),  handles.img_size(3));
    
    for c_roi = 1:numel(active_rois) %use meshgrid if this is slow
        h_roi = active_rois(c_roi);
        oef_slice = active_roi_slices(c_roi);
        for x = 1:img_size(1)
            for y = 1:img_size(2)
                if inROI(h_roi, x,y)
                    oef_mask(y,x,oef_slice) = true;
                end
            end
        end
    end
    
    guidata(hObject, handles)
    
    wrap_text = 'Calibrating OEF calculation. Please wait...';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    drawnow
    figure(handles.p_bar)
    waitbar(6/7, handles.p_bar, 'Calibrating OEF calculation: user input accepted...')
    % set(handles.figure1, 'Pointer', 'watch')
    
    %% calibrate k for OEF WM = OEF_WM
    k_upper = 10000000;
    k_lower = 0;
    incomplete = true;
    %open log file
    fileID = fopen(strcat(target_folder, '/OEF_calibration_log.txt'), 'a');
    fprintf(fileID, "\n");
    fprintf(fileID, "Begin logging OEF calibration: \n");
    
    fail_counter = 0;
    do_BzD = handles.BzD;
    do_SVD = handles.do_SVD;
    previous_tmp_k = -1;  
    while incomplete
        %tau is the capillary transit time
        
        % k       = 118;      % [1/s]
        
        tmp_k = k_lower + round((k_upper-k_lower)/2, 1);
        
        if tmp_k == previous_tmp_k
            fprintf(fileID, "OEF Calibration failed. Required OEF can not be achieved with k in the range [0 1E7]" + "\n");
            fprintf(fileID, "Proceeding with a default value of k = 50 /s" + "\n");
            message = "OEF Calibration failed. Required OEF can not be achieved with k in the range [0 1E7]";
            disp(message);
            opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
            uiwait(errordlg('\fontsize{10}OEF calibration failed. Required OEF can not be achieved with k in the range [0 1E7]', 'DECONVOLER: OEF calibration', opts));
            tmp_k = 50;
            break;
        end
        
        B       = 0.1943;   % ml/ml
        c_a     = 0.95*B;   % ml/ml
        c_0     = c_a;
        
        xspan   = [0 1];
        taus    = (0:0.1:100)';
        dtau    = taus(2)-taus(1);
        
        options = odeset('AbsTol',1e-9, 'RelTol', 1e-4); %% set solver options
        Q = zeros(length(taus),1);
        
        for i = 1:length(taus)
            
            tau = taus(i);
            
            odefun = @(x, C)dc_dx(x,C, tau, tmp_k, Pt, P50, h);
            [~,C] = ode23tb(odefun, xspan,c_0,options); %% solve equations
            %@@@ integrates the system of differential equations C' = f(C,x) from
            %tspan(0) to tspan(1), with initial conditions given by c_0
            Q(i) = 1 - C(end)/C(1);
        end
        
        % save("Q_right.mat", "Q", "taus")
        %
        oef_BzD   = NaN(img_size(1:3));
        oef_svd   = NaN(img_size(1:3));
        xrange = slice_range(1):slice_range(2);
        yrange = slice_range(3):slice_range(4);
        zrange = slice_range(5):slice_range(6);
        for x = xrange
            tmp_oef_BzD = NaN(img_size(2:3));
            tmp_oef_svd = NaN(img_size(2:3));
            for y = yrange
                for z = zrange
                    if  oef_mask(x,y,z)
                        if do_BzD
                            r = bezier_residue_function(fitd_omega(x,y,z,:), taus);
                            H = -diff(r)/dtau;
                            H(end+1) = 0;
                            %                             h = -customDiff(r)/dtau;
                            tmp_oef_BzD(y,z) = 100*trapz(taus,H.*Q);
                        end
                        if do_SVD
                            try
                                r = interp1(t, fitd_r_svd(x,y,z,:), taus);
                            catch
                                r = zeros(size(taus));
                            end
                            %                             h = -customDiff(r)/dtau;
                            H = -diff(r)/dtau;
                            H(end+1) = 0;
                            tmp_oef_svd(y,z) = 100*trapz(taus,H.*Q);
                        end
                    end
                end
            end
            oef_BzD(x,:,:) = tmp_oef_BzD;
            oef_svd(x,:,:) = tmp_oef_svd;
        end
        
        
        
        if do_SVD; oef_svd(isnan(oef_svd)) = []; wm_oef = mean(oef_svd); end
        if do_BzD; oef_BzD(isnan(oef_BzD)) = []; wm_oef = mean(oef_BzD); end
        
        if isfinite(wm_oef)&&(wm_oef>0)&&isreal(wm_oef)
            if (wm_oef) > 1.05*OEF_WM
                k_upper = tmp_k;
                fprintf(fileID, "k = " + num2str(tmp_k) + " is too high. OEF WM = " + num2str(wm_oef) + "\n");
                disp("k = " + num2str(tmp_k) + " is too high. OEF = " + num2str(wm_oef));
            elseif (wm_oef) < OEF_WM
                k_lower = tmp_k;
                fprintf(fileID, "k = " + num2str(tmp_k) + " is too low. OEF WM = " + num2str(wm_oef) + "\n");
                disp("k = " + num2str(tmp_k) + " is too low. OEF = " + num2str(wm_oef));
            else
                incomplete = false;
                fprintf(fileID, "k = " + num2str(tmp_k) + " is perfect.  OEF WM = " + num2str(wm_oef) + "\n");
                disp("k = " + num2str(tmp_k) + " is perfect.  OEF = " + num2str(wm_oef));
            end
        else
            
            fprintf(fileID, "k = " + num2str(tmp_k) + " : OEF WM is not finite" + "\n");
            k_upper = 0.5*k_upper;
            fail_counter = fail_counter+1;
            if fail_counter > 20
                fprintf(fileID, "OEF Calibration failed. Proceeding with default value of k = 50" + "\n");
                 message = "OEF Calibration failed. Proceeding with default value of k = 50";
                 disp(message);
            opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
            uiwait(errordlg('\fontsize{10}OEF calibration failed. Calibration settings gave undefined OEF.', 'DECONVOLER: OEF calibration', opts));
                tmp_k = 50;
                break
            end
            disp("k = " + num2str(tmp_k) + " : OEF not finite")
        end
        
        previous_tmp_k = tmp_k; %use this to identify when required k is outside maximum range        
    end
    
    fprintf(fileID, "End logging OEF calibration. \n");
    fclose(fileID);
    k = tmp_k;
    
    oef_rois = active_rois;
    oef_roi_slices = active_roi_slices;
    save(fullfile(handles.target_folder, "OEF_calibration.mat"), "oef_rois", "oef_roi_slices", "k", "Pt", "OEF_WM", "P50", "h")
    tmp = get(handles.edit14, 'String');
    tmp{end+1} = '- - - - OEF calibration file - - - -'; tmp{end+1} = fullfile(handles.target_folder, 'OEF_calibration.mat');
    set(handles.edit14, 'String',tmp)
    guidata(hObject, handles)
end

%% before computing OEF, ask user to verify or change k, Pt, P50 and h
%ask user for Pt here again (Pt is tissue oxygen tension in mmHg),  P50
%(P50 is oxygen pressure required for half saturation) and h (Hill coefficient) 
wrap_text = 'Waiting for user input...';
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
opts.Interpreter = 'tex';
prompt = {'\color{blue} \fontsize{10} Calibrated rate constant, k [s^{-1}] :', '\color{blue} \fontsize{10} Tissue oxygen tension, Pt [mmHg] :'...
    , '\color{blue} \fontsize{10} Oxygen pressure required for half saturation, P50 [mmHg] :', '\color{blue} \fontsize{10} Hill coefficient, h :'};
dlgtitle = 'DECONVOLVER: Confirm/edit OEF calculation settings';
dims = [1 90];
defaults = {num2str(k), num2str(Pt), num2str(P50), num2str(h)};
response = inputdlg(prompt, dlgtitle, dims, defaults, opts);
if isempty(response) %if use clicks cancel, do nothing
    wrap_text = 'Proceeding with calibrated values...';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
else
    tmp_k = str2double(response{1});
    tmp_Pt = str2double(response{2});
    tmp_P50 = str2double(response{3});
    tmp_h = str2double(response{4});
    if isfinite(tmp_k)&&isfinite(tmp_Pt)&&isfinite(tmp_P50)&&isfinite(tmp_h)
        wrap_text = 'User-supplied OEF calculation settings saved.';
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
        k = tmp_k;
        Pt = tmp_Pt;
        P50 = tmp_P50;
        h = tmp_h;
    else
        wrap_text = 'Invalid input. Proceeding with calibrated values...';
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    end
end


%% write final values of k, Pt and P50 to file before continuing
fileID = fopen(strcat(target_folder, '/OEF_calibration_log.txt'), 'a');
fprintf(fileID, "\n");
fprintf(fileID, "OEF calculation settings: \n");
fprintf(fileID, "k = " + num2str(k) + "/s , Pt = " + num2str(Pt) +  " mmHg, P50 = " + num2str(P50) + " mmHg and h = " + num2str(h) + "\n");
fclose(fileID);

%% now compute actual OEF
report('Computing OEF...', handles)
figure(handles.p_bar)
waitbar(6/7, handles.p_bar, 'Calculating OEF and CMRO2...')
% set(handles.figure1, 'Pointer', 'watch')
mask = handles.mask;
%tau is the capillary transit time

B       = 0.1943;   % ml/ml
c_a     = 0.95*B;   % ml/ml
c_0     = c_a;

xspan   =  [0 1]; %normalised distance
taus    = (0:0.1:100)';
dtau    = taus(2)-taus(1);

options = odeset('AbsTol',1e-9, 'RelTol', 1e-4); %% set solver options
Q = zeros(length(taus),1);

for i = 1:length(taus)
    
    tau = taus(i);
    
    odefun = @(x, C)dc_dx(x,C, tau, k, Pt, P50, h);
    [~,C] = ode23tb(odefun, xspan,c_0,options); %% solve equations
    %@@@ integrates the system of differential equations C' = f(C,x) from
    %tspan(0) to tspan(1), with initial conditions given by c_0
    Q(i) = 1 - C(end)/C(1);
end
% save("Q_right.mat", "Q", "taus")
oef_BzD   = zeros(img_size(1:3));
cmro2_BzD = zeros(img_size(1:3));
oef_svd   = zeros(img_size(1:3));
cmro2_svd = zeros(img_size(1:3));

xrange = slice_range(1):slice_range(2);
yrange = slice_range(3):slice_range(4);
zrange = slice_range(5):slice_range(6);
do_BzD = handles.BzD;
do_SVD = handles.do_SVD;

for x = xrange
    tmp_oef_BzD = zeros(img_size(2:3));
    tmp_cmro2_BzD = zeros(img_size(2:3));
    tmp_oef_svd = zeros(img_size(2:3));
    tmp_cmro2_svd = zeros(img_size(2:3));
    
    for y = yrange
        for z = zrange
            if mask(x,y,z)
                if do_BzD
                    r = bezier_residue_function(fitd_omega(x,y,z,:), taus);
                    H = -diff(r)/dtau;
                    H(end+1) = 0;
                    tmp_oef_BzD(y,z) = 100*trapz(taus,H.*Q);
                    tmp_cmro2_BzD(y,z) = c_a * tmp_oef_BzD(y,z) * fitd_cbf(x,y,z);
                end
                if do_SVD
                    try
                        r = interp1(t, fitd_r_svd(x,y,z,:), taus);
                    catch
                        r = zeros(size(taus));
                    end
                    %                     h = -customDiff(r)/dtau;
                    H = -diff(r)/dtau;
                    H(end+1) = 0;
                    tmp_oef_svd(y,z) = 100*trapz(taus(1:numel(H)),H.*Q);
                    tmp_cmro2_svd(y,z) = c_a * tmp_oef_svd(y,z) * fitd_cbf_svd(x,y,z);
                end
            end
        end
    end
    oef_BzD(x,:,:)   = tmp_oef_BzD;
    cmro2_BzD(x,:,:)  = tmp_cmro2_BzD;
    oef_svd(x,:,:)   = tmp_oef_svd;
    cmro2_svd(x,:,:) = tmp_cmro2_svd;
end



%save results
try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
save_data.this_folder = target_folder; save_data.this_format = handles.save_format;
if handles.BzD; save_data.this_name = 'oef_BzD'; save_data.data_to_save = oef_BzD; save_this_file(save_data); end
if handles.do_SVD; save_data.this_name = 'oef_SVD'; save_data.data_to_save = oef_svd; save_this_file(save_data); end
if handles.BzD; save_data.this_name = 'cmro2_BzD'; save_data.data_to_save = cmro2_BzD; save_this_file(save_data); end
if handles.do_SVD; save_data.this_name = 'cmro2_SVD'; save_data.data_to_save = cmro2_svd; save_this_file(save_data); end

% oef_right = squeeze(oef_BzD(:,:,9));
% save("oef_right_k" + num2str(k)+".mat", "oef_right")
% oef_wrong = squeeze(oef_BzD(:,:,9));
% save("oef_wrong_k" + num2str(k)+".mat", "oef_wrong")

report('OEF calculation complete', handles)
set(handles.figure1, 'Pointer', 'arrow')
set(handles.slider1, 'Callback', {@slider1_Callback, hObject})
set(handles.slider8, 'Callback', {@slider8_Callback, hObject})
guidata(hObject, handles)
end


%% Helper functions
function prepare_sliders(hObject, ~, ~)
handles = guidata(hObject);
initial_slice = round(handles.img_size(3)/2);
t_min_signal = handles.t_min_signal;

initial_slice_image = squeeze(handles.dsc_data_filtered(:,:,initial_slice,handles.t_min_signal));
clim_im = [0   max(initial_slice_image(isfinite(initial_slice_image)))]; %color limit, can be changed using sliders
handles.dsc_data_image = imshow(zeros(handles.img_size(1), handles.img_size(2)),clim_im, 'colormap', gray, 'Parent', handles.axes1);
title(handles.axes1, ['Slice ' num2str(initial_slice)])
set(handles.dsc_data_image, 'CData', initial_slice_image ) %display initial slice
set(handles.slider1, 'Min', 1, 'Max', handles.img_size(3), ... %setup slider for changing slice
    'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', initial_slice)
set(handles.slider8, 'Visible', 'on')%horizontal scrolling on axes1
set(handles.slider1, 'Visible', 'on')%horizontal scrolling on axes1
set(handles.slider8, 'Min', 1, 'Max', handles.img_size(4), ...
    'SliderStep', [1, 1]/(handles.img_size(4) - 1), 'Value', t_min_signal)
clim_max = double(max(handles.dsc_data_filtered(isfinite(handles.dsc_data_filtered)))) + 1E-10;
if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/999; end
if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end %setup slider for changing window setting
set(handles.slider6, 'Visible', 'on') %color setting
set(handles.slider6, 'Min', 10E-10, 'Max', clim_max, ...
    'SliderStep', [step_small, step_big], 'Value', clim_im(2))
set(handles.slider10, 'Visible', 'on') %color setting
set(handles.slider10, 'Min', 0, 'Max', clim_max-1E-10, ...
    'SliderStep', [step_small, step_big], 'Value', 0)
% hObject = handles.output;
guidata(hObject, handles)
end

function initialise_roi_placement(hObject, ~)
handles = guidata(hObject);
max_N_roi = 128*128;
handles.roi_counter = 0;
handles.rois = gobjects(max_N_roi,1);
handles.roi_slices = zeros(max_N_roi, 1);
handles.roi_mask = false(handles.img_size(1), handles.img_size(2), handles.img_size(3));
guidata(hObject, handles)
end

function Proceed(hObject, ~, ~)
handles = guidata(hObject);
set(handles.pushbutton26, 'Visible', 'off');
set(handles.pushbutton24, 'Visible', 'off');
set(handles.pushbutton24, 'Callback', [])
set(handles.pushbutton26, 'Callback', [])
set(handles.pushbutton30, 'String', 'Exit AIF selection')
set(handles.pushbutton30, 'Callback', [])
set(handles.slider1, 'Callback', []);
set(handles.slider8, 'Callback', []);
set(handles.pushbutton30, 'Visible', 'off')
guidata(hObject, handles)
end

function addROI(hObject, ~, ~)
handles = guidata(hObject);
roi = drawfreehand(handles.axes1, 'Color', 'r', 'LineWidth',1, 'FaceAlpha', 0);
handles.roi_counter = handles.roi_counter+1;
handles.rois(handles.roi_counter) = roi;
handles.roi_slices(handles.roi_counter) = get(handles.slider1, 'value');
guidata(hObject, handles) %update master copy
end

function change_slice(hObject, ~, ~)
handles = guidata(hObject);
current_slice = round(get(hObject, 'Value'));
current_time_point = round(get(handles.slider8, 'Value'));
set(handles.dsc_data_image, 'CData', squeeze(handles.dsc_data_filtered(:,:,current_slice,current_time_point)))
active_rois = handles.rois(ishandle(handles.rois));
active_roi_slices = handles.roi_slices(ishandle(handles.rois));
set(active_rois, 'Parent', []);
set(active_rois(active_roi_slices == current_slice), 'Parent', handles.axes1);
% title(handles.axes1, ['Signal: Slice ' num2str(current_slice) ' : Time point ' num2str(current_time_point)])
title(handles.axes1, ['Slice ' num2str(current_slice)])
guidata(hObject, handles)
end

function change_timepoint(hObject, ~ , ~)
handles = guidata(hObject);
time_point = round(get(hObject,'Value'));
current_slice = round(get(handles.slider1,'Value'));
set(handles.dsc_data_image, 'CData', squeeze(handles.dsc_data_filtered(:,:,current_slice, time_point))  )
% title(handles.axes1, ['Signal' ': Slice ' num2str(current_slice) ' : Time-point ' num2str(time_point)],...
%     'Interpreter', 'None')
guidata(hObject, handles)
end

function prepare_dsc_data(hObject)
handles = guidata(hObject);
if isfield(handles, 'mask') %if mask has been generated
    handles.dsc_data_filtered = handles.dsc_data.*0;
    for sl = 1:handles.img_size(3) %filter image to be displayed using the mask
        for tp = 1:handles.img_size(4)
            handles.dsc_data_filtered(:,:,sl, tp) = squeeze(handles.dsc_data(:,:,sl, tp)).*squeeze(handles.mask(:,:,sl));
        end
    end
else
    handles.dsc_data_filtered = handles.dsc_data;
end
guidata(hObject, handles)
end

function clearROIs(hObject, ~, ~)
%clear all drawn ROIS
handles = guidata(hObject);
clim_im = get(handles.axes1, 'CLim');
active_rois = handles.rois(ishandle(handles.rois));
% active_roi_slices = handles.roi_slices(ishandle(handles.rois));
for c_ar = 1:numel(active_rois)
    delete(active_rois(c_ar))
end
current_slice = round(get(handles.slider1, 'Value'));
current_timepoint = round(get(handles.slider8, 'Value'));
cla(handles.axes1, 'reset')
handles.dsc_data_image = imshow(zeros(handles.img_size(1), handles.img_size(2)), clim_im, 'colormap', gray, 'Parent', handles.axes1);
set(handles.dsc_data_image, 'CData', (squeeze(handles.dsc_data_filtered(:,:,current_slice,current_timepoint))) )
guidata(hObject, handles)
initialise_roi_placement(hObject)
end

