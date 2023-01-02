function varargout = DECONVOLVER(varargin)
%      DECONVOLVER performs analysis of DSC-MRI data using either
%      deconvolution with Bezier curves (BzD) or SVD-based methods. BzD is
%      available with and without correction for delay and/or dispersion
%      effects. All variants of SVD deconvolution are availabe: sSVD, cSVD
%      and oSVD.
%
%      The input data to DECONVOLVER is a 4-dimensional DSC time-series,
%      preferably in the format NIFTI (*.nii or *.nii.gz). The capability
%      of the program to load a 4-D volume from a series of Dicom images is
%      still under development.
%
%      DECONVOLVER is equipped with AIF and VOF selection algorithms.
%      There are three variants of these algorithms: Automatic,
%      Semi-automatic and Manual.
%
%      DECONVOLVER computes the quantities CBF, CBV, MTT, OEF, CMRO2, TTP,
%      R10, R50 and Residue functions and saves them as per user
%      specification.
%
%      Arbitrary data can also be loaded onto the axes and saved from them
%      independent of any analysis.
%
%      For more information and help, contact the author:
%      Arthur Chakwizira
%      arthur.chakwizira@med.lu.se
%      Medical Radiation Physics, Lund University, Sweden
%
%      Part of the code base is an extension of preliminary work by:
%      Andre Ahlgren

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @DECONVOLVER_OpeningFcn, ...
    'gui_OutputFcn',  @DECONVOLVER_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DECONVOLVER is made visible.
function DECONVOLVER_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DECONVOLVER (see VARARGIN)
%--------------------------------------------------------------------------

%Enable logging
clc
log_fn = 'DECONVOLVER_log.txt';
fileID = fopen(log_fn, 'w');
if fileID ~= -1; fclose(fileID); diary(log_fn); end


% DECONVOLVER needs the functions in folder named 'functions' to run
% path_to_functions = fullfile(pwd, 'functions');
% if ~isfolder(path_to_functions) %if functions folder is non-existent
%     opts.Interpreter = 'tex'; opts.WindowStyle = 'modal'; %raise an error and terminate execution
%     errordlg('\fontsize{10} \color{red}All functions needed for DECONVOLVER to run are missing! DECONVOLVER can not run.', 'Fatal error', opts);
%     error('FATAL ERROR! Functions needed for DECONVOLVER to run are missing!');
% end
% addpath(path_to_functions) %add path to functions


%Setting default values for the program.
movegui(gcf, 'center') %move the window to the centre of the screen when it opens
handles.axes1_vaccant = false; %some functions need to know which axes is not currently displaying data
handles.axes2_vaccant = false;
handles.cmap_index = 0; %initialise colormap selection index
handles.colormaps = {'hot', 'jet', 'parula', 'hsv', 'cool', 'gray'}; %possible colormaps
set(handles.slider1, 'Visible', 'off') %slice change on axes1
set(handles.slider2, 'Visible', 'off'); %slice change on axes2
set(handles.slider6, 'Visible', 'off') %upper window setting on axes1
set(handles.slider10, 'Visible', 'off')%lower window setting on axes1
set(handles.slider7, 'Visible', 'off')%upper window setting on axes2
set(handles.slider11, 'Visible', 'off')%lower window setting on axes2
set(handles.slider8, 'Visible', 'off')%horizontal scrolling on axes1
set(handles.slider9, 'Visible', 'off')%horizontal scrolling on axes2
set(handles.pushbutton15, 'Visible', 'off')%colormap on axes1
set(handles.pushbutton16, 'Visible', 'off')%colormap on axes2
set(handles.togglebutton2, 'Visible', 'off')%colorbar on axes1
set(handles.togglebutton3, 'Visible', 'off')%colorbar on axes2
set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes2
handles.data_onto_axes1 = false; %is used to alert objects that user has loaded data onto axes1
handles.data_onto_axes2 = false; %is used to alert objects that user has loaded data onto axes2
handles.plotting_4d_data = false; %alert objects that user is visualising 4D data as plots
handles.ext_plot = 0; %keeps track of external figures opened to allow saving of axes data in formats png,jpg,etc.
handles.calculate_oef = false; %calculate_oef
handles.calculate_cmro2 = false; %calculate_cmro2
handles.calculate_r10 = false; %calculate_r10
handles.calculate_r50 = false; %calculate_r50
handles.data_is_for_aif = false; %alert objects that user is currently selecting AIF
handles.data_is_for_vof = false; %alert objects that user is currently selecting VOF
%all the variables below are used to check if some function has reserved a given slider
handles.cbf_slider = 0;
handles.mtt_slider = 0;
%_________________________________________________________________________
set(handles.axes1,'XTick', [], 'YTick', []) %set default appearance for axes
set(handles.axes2,'XTick', [], 'YTick', [])
set(handles.axes1, 'box', 'on')
set(handles.axes2, 'box', 'on')
%__________________________________________________________________________
handles.kappa = 0.705; %kappa; the hematocrit constant
handles.calculate_aif = false; %will invoke AIF selection algorithm if true
handles.calculate_vof = false; %will invoke AIF selection algorithm if true
handles.include_vof = false; %choose whether to include a VOF for PVE correction
handles.BzD = true; handles.with_delay = true; handles.with_dispersion = false; %initial deconvolution algorithm settings
handles.do_SVD = false; handles.do_oSVD = false; handles.do_cSVD = false; handles.do_sSVD = false; %deconvolution algorithm settings
set(handles.checkbox26, 'value', 1) %calculate_oef
set(handles.checkbox28, 'value', 0) %calculate_r10
%the following are initial settings for what is displayed on the axes; Default is CBF on left axes and MTT on right axes)
handles.show_cbf = true; handles.cbf_axes = handles.axes1; handles.cbf_slider = handles.slider1; handles.cbf_window_slider_up = handles.slider6; handles.cbf_window_slider_low = handles.slider10;
handles.show_mtt = true; handles.mtt_axes = handles.axes2; handles.mtt_slider = handles.slider2; handles.mtt_window_slider_up = handles.slider7; handles.mtt_window_slider_low = handles.slider11;
%__________________________________________________________________________
handles.save_format = 'NIFTI (*.nii)'; %save_format, default is NIFTI
wrap_text = 'WELCOME'; %default text shown on display panel
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
%v 1.4
handles.dsc_data_loaded = false; %true when data has been loaded to working memory
handles.aif_exists = false;
%Select file filter to be used for saving of files (needed for AIF / VOF)
file_filter = {'NIFTI *.nii'; 'DICOM *.dcm' ; 'MATLAB *.mat'; 'Text *.txt'; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat' };
handles.file_filter = file_filter;
handles.mask_threshold = 100;
handles.hObject = handles.figure1;
set([handles.pushbutton24, handles.pushbutton25, handles.pushbutton26, handles.pushbutton27, handles.pushbutton30], 'Visible', 'off')
set([handles.pushbutton31, handles.pushbutton32, handles.pushbutton33, handles.pushbutton34, handles.pushbutton35], 'Visible', 'off') %same buttons on axes2

% handles.save_format = %TODO
% Choose default command line output for DECONVOLVER
handles.output = handles.figure1;
handles.hObject = handles.figure1;
 handles.save_results = true;
 handles.notify_when_done= false;
handles.slider1_Callback = handles.slider1.Callback;
handles.slider8_Callback = handles.slider8.Callback;
handles.choosing_region_to_analyse = false;
% set(handles.pushbutton23, 'Enable', 'off')
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = DECONVOLVER_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in AIF SELECTION METHOD.
function popupmenu4_Callback(hObject, eventdata, handles)
%Manages the user's choice of AIF selection algorithm
aif_choices = cellstr(get(hObject,'String'));
aif_choice = aif_choices{get(hObject,'Value')};
switch aif_choice
    case 'Load'
        handles.calculate_aif = false;
    case 'Select: Semi-Automatic'
        handles.calculate_aif = true;
        handles.calculate_aif_semi_auto = true;
        handles.calculate_aif_auto = false;
        handles.calculate_aif_manual = false;
    case 'Select: Automatic'
        handles.calculate_aif = true;
        handles.calculate_aif_auto = true;
        handles.calculate_aif_semi_auto = false;
        handles.calculate_aif_manual = false;
    case 'Select: Manual'
        handles.calculate_aif = true;
        handles.calculate_aif_auto = false;
        handles.calculate_aif_semi_auto = false;
        handles.calculate_aif_manual = true;
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in VOF SELECTION METHOD.
function popupmenu5_Callback(hObject, eventdata, handles)
%Manages the user's choice of VOF selection algorithm
vof_choices = cellstr(get(hObject,'String'));
vof_choice = vof_choices{get(hObject,'Value')};
switch vof_choice
    case 'None'
        handles.include_vof = false;
    case 'Load'
        handles.include_vof = true;
        handles.calculate_vof = false;
    case 'Select: Semi-Automatic'
        handles.include_vof = true;
        handles.calculate_vof = true;
        handles.calculate_vof_semi_auto = true;
        handles.calculate_vof_auto = false;
        handles.calculate_vof_manual = false;
    case 'Select: Automatic'
        handles.include_vof = true;
        handles.calculate_vof = true;
        handles.calculate_vof_auto = true;
        handles.calculate_vof_semi_auto = false;
        handles.calculate_vof_manual = false;
    case 'Select: Manual'
        handles.include_vof = true;
        handles.calculate_vof = true;
        handles.calculate_vof_auto = false;
        handles.calculate_vof_semi_auto = false;
        handles.calculate_vof_manual = true;
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when user presses RUN BUTTON.
function pushbutton4_Callback(hObject, eventdata, handles)
%THIS IS THE MAIN PART OF THE PROGRAM

%% ANALYSIS OF DSC-MRI DATA
hObject = handles.hObject;
handles = guidata(hObject);
%% CHECK IF AIF and DATA EXIST
if ~handles.dsc_data_loaded;  wrap_text = 'No DSC-MRI data loaded. Load data to proceed.';  set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r'); return; end
if ~handles.aif_exists;  wrap_text = 'AIF not selected. Load or select AIF to proceed.';  set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r'); return; end


%% GET CONFIGURATION SETTINGS <>
handles.data_is_for_aif = false;%if true, slider1 will assume that concentration data is being shown on axes1
%restore window to its state when it was opened
set(handles.slider1, 'Visible', 'off') %slice change on axes1
set(handles.slider2, 'Visible', 'off'); %slice change on axes2
set(handles.slider6, 'Visible', 'off') %upper window setting on axes1
set(handles.slider10, 'Visible', 'off')%lower window setting on axes1
set(handles.slider7, 'Visible', 'off')%upper window setting on axes2
set(handles.slider11, 'Visible', 'off')%lower window setting on axes2
set(handles.slider8, 'Visible', 'off')%horizontal scrolling on axes1
set(handles.slider9, 'Visible', 'off')%horizontal scrolling on axes2
set(handles.pushbutton15, 'Visible', 'off')%colormap on axes1
set(handles.pushbutton16, 'Visible', 'off')%colormap on axes2
set(handles.togglebutton2, 'Visible', 'off')%colorbar on axes1
set(handles.togglebutton3, 'Visible', 'off')%colorbar on axes2
set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes2
%__________________________________________________________________________
set(handles.edit7, 'String', '') %clear display/message panel
handles.data_onto_axes1 = false; %alerts objects if user has loaded data onto axes
handles.data_onto_axes2 = false;
handles.data_is_for_aif = false; %alerts objects if user is currently selecting AIF
handles.data_is_for_vof = false; %alerts objects if user is currently selecting VOF
handles.calculate_oef = get(handles.checkbox26, 'Value'); %calculate_oef
handles.calculate_r10 = get(handles.checkbox28, 'Value'); %calculate_r10
%restore axes to their state when window opened
cla(handles.axes1, 'reset')
cla(handles.axes2, 'reset')
set(handles.axes1,'XTick', [], 'YTick', [])
set(handles.axes2,'XTick', [], 'YTick', [])
set(handles.axes1, 'box', 'on')
set(handles.axes2, 'box', 'on')
set([handles.pushbutton31, handles.pushbutton32,handles.pushbutton33, handles.pushbutton34, handles.pushbutton35], 'Visible', 'off') %initially turn off all these buttons
set([handles.pushbutton24, handles.pushbutton25,handles.pushbutton26, handles.pushbutton27, handles.pushbutton30], 'Visible', 'off') %initially turn off all these buttons
reset_axes(handles)

if isfield(handles, 'old_mask'); handles.mask = handles.old_mask; end %undo restrict analysis to user roi
guidata(hObject, handles) %update master

%% GET VOF
if handles.include_vof %if user wants to include VOF (for PVE correction)
   if ~handles.vof_exists;  wrap_text = 'VOF not selected. Select or load VOF to proceed.';  set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r'); return; end 
end

%% RESTRICT ANALYSIS TO SPECIFIC REGIONS?
%Invoke question dialog to ask user if user wants to define ROI to restrict
%analysis to particular regions of the brain.
opts.Interpreter = 'tex';
opts.Default = 'No';
response = questdlg('\fontsize{10}Program allows ROI placement to restrict analysis to user-defined regions. Continue with ROI placement?', ...
    'DECONVOLVER', 'No', 'Yes','Stop analysis', opts);
switch response
    case 'No'; wants_roi = false;
    case 'Yes'; wants_roi = true;
    case 'Stop analysis'; set(handles.edit7, 'String', 'Analysis terminated by user.', 'ForeGroundColor', 'r'); return
    case ''; set(handles.edit7, 'String', 'Analysis terminated by user.', 'ForeGroundColor', 'r'); return
end

if wants_roi; restrict_analysis_to_user_defined_regions(hObject); end
handles = guidata(hObject); %fetch updated handles

%generate results folder
handles.target_folder = generate_target_folder(handles);
guidata(hObject, handles) %update master
%% BEZIER CURVE DECONVOLUTION
pause(0.5)
if handles.BzD %if user chose Bezier Curve deconvolution
        handles = get_deconvolution_settings(handles, 'bzd'); %get these settings
        if handles.user_cancelled; handles.user_cancelled = []; return; end %if user cancelled input of these settings, terminate execution
        guidata(hObject, handles) %update app data
    
    %     handles.p_bar = waitbar(0, ['Bezier-curve deconvolution: ' '0% complete'], 'Name', 'DECONVOLVER'); %create progress bar
    
    wrap_text = 'Running Bezier-curve deconvolution ...'; %display this message on the main window
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    pause(0.5)
    childlist = get(handles.figure1, 'Children');
    for c = 1:numel(childlist)
        try set(childlist(c), 'Enable', 'off'); catch; end
    end
    %         set(handles.figure1, 'pointer', 'watch')
    pause(0.5)
    params = bezier_deconvolve_gui(handles); %otherwise without live data visualisation
    fitd_omega     = params(:,:,:,1:5); %control points
    fitd_cbf      	= params(:,:,:,6).*6000; %CBF in ml/100g/min
    %get delay and dispersion if relevant
    %__________________________________________________________________________
    if handles.with_dispersion && ~handles.with_delay;   fitd_dk  =  params(:,:,:,7:8);  fitd_delay = false;  end
    if handles.with_delay && ~handles.with_dispersion;   fitd_delay = params(:,:,:,7);  fitd_dk = false;  end
    if handles.with_delay && handles.with_dispersion
        fitd_delay = params(:,:,:,7);
        fitd_dk  =  params(:,:,:,8:9);
    end
    if ~handles.with_delay && ~handles.with_dispersion
        fitd_delay  	= false;
        fitd_dk = false;
    end
    %__________________________________________________________________________
    if handles.include_vof %if user chose to include VOF for PVE correction
        fitd_cbf = fitd_cbf*handles.pve_correction; %incorperate correction
    end
    %save results______________________________________________________________
    try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
    save_data.this_folder = handles.target_folder; save_data.this_format = handles.save_format;
    save_data.this_name = 'control_points_BzD'; save_data.data_to_save = fitd_omega; save_this_file(save_data)
    save_data.this_name = 'cbf_Bzd'; save_data.data_to_save = fitd_cbf; save_this_file(save_data)
    if handles.with_delay; save_data.this_name = 'delay_BzD'; save_data.data_to_save = fitd_delay; save_this_file(save_data); end
    if handles.with_dispersion; save_data.this_name = 'disp_kernel_BzD'; save_data.data_to_save = fitd_dk; save_this_file(save_data); end
else
    fitd_cbf = false; fitd_delay = false; fitd_dk = false; fitd_omega = false;
end
%__________________________________________________________________________
%% SVD deconvolution
if handles.do_SVD %if user chose one of the SVD methods
        if handles.do_oSVD; caller = 'osvd'; end %let 'get_deconvolution_settings' which algorithm has been chosen
        if handles.do_cSVD; caller = 'csvd'; end
        if handles.do_sSVD; caller = 'ssvd'; end
        handles = get_deconvolution_settings(handles, caller); %get the relevant settings
        if handles.user_cancelled; handles.user_cancelled = false; return; end %if user cancelled input dialog, return
        guidata(hObject, handles)
    pause(0.5)
    %create progress bar with correct text depending on chosen algorithm____
    if handles.do_oSVD; wrap_text = 'Running oSVD deconvolution ...'; end
    if handles.do_cSVD; wrap_text = 'Running cSVD deconvolution ...'; end
    if handles.do_sSVD; wrap_text = 'Running sSVD deconvolution ...';  end
    %_________________________________________________________________________
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b') %display same message on main window
    pause(0.5)
    %     handles.p_bar = waitbar(0, [p_bar_id '0% complete'], 'Name', 'DECONVOLVER');
    
%     childlist = get(handles.figure1, 'Children');
%     for c = 1:numel(childlist)
%         try set(childlist(c), 'Enable', 'off'); catch; end
%     end
    %         set(handles.figure1, 'pointer', 'watch')
    
    pause(0.5)

        if handles.do_oSVD
            [fitd_r_svd, fitd_cbf_svd, fitd_delay_svd] = oSVD_deconvolution_gui(handles); %otherwise call the conventional function
        end
        if handles.do_cSVD
            [fitd_r_svd, fitd_cbf_svd, fitd_delay_svd] = cSVD_deconvolution_gui(handles); %otherwise call the conventional function
        end
        if handles.do_sSVD
            [fitd_r_svd, fitd_cbf_svd, fitd_delay_svd] = sSVD_deconvolution_gui(handles); %otherwise call the conventional function
        end
    %__________________________________________________________________________
    if handles.include_vof %include vof for PVE correction
        fitd_cbf_svd = fitd_cbf_svd*handles.pve_correction;
    end
    %save results______________________________________________________________
    try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
    save_data.this_folder = handles.target_folder; save_data.this_format = handles.save_format;
    save_data.this_name = 'cbf_SVD'; save_data.data_to_save = fitd_cbf_svd; save_this_file(save_data)
else
    fitd_cbf_svd = false; fitd_r_svd = false; fitd_delay_svd = false;
end

%save results to handles structure_________________________________________
handles.fitd_cbf = fitd_cbf;
handles.fitd_cbf_svd = fitd_cbf_svd;
handles.fitd_omega = fitd_omega;
handles.fitd_r_svd = fitd_r_svd;
handles.fitd_dk = fitd_dk;
handles.fitd_delay = fitd_delay;
handles.fitd_delay_svd = fitd_delay_svd;
%__________________________________________________________________________
%% CALCULATE CBV and MTT, AS WELL AS R10/R50, OEF/CMRO2 AND TTP
%CBV, deconvolution-independent____________________________________________
set(handles.edit7, 'String', 'Computing additional quantities...' ,'ForeGroundColor', 'blue')
handles.p_bar = waitbar(0, 'Calculating CBV: {CBV = \kappa_H \cdot \int{C(t)}/ \int{AIF} }', 'Name', 'DECONVOLVER'); %progress bar
cbv_without_mtt = cbv_no_mtt_gui(handles); %ml/100g
if handles.include_vof; cbv_without_mtt = cbv_without_mtt.*handles.pve_correction; end
try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
save_data.this_folder = handles.target_folder; save_data.this_format = handles.save_format;
save_data.this_name = 'cbv_without_deconvolution'; save_data.data_to_save = cbv_without_mtt; save_this_file(save_data)
handles.cbv_without_mtt = cbv_without_mtt; %save to handles

%Bezier residue functions__________________________________________________
waitbar(1/7, handles.p_bar, 'Converting Bezier control points to residue functions')
if handles.BzD; r_BzD = get_r_from_cp(handles);
    try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
    save_data.this_folder = handles.target_folder; save_data.this_format = handles.save_format;
    save_data.this_name = 'residue_functions_BzD'; save_data.data_to_save = r_BzD; save_this_file(save_data)
end

%MTT as CBV/CBF____________________________________________________________
waitbar(2/7, handles.p_bar, 'Calculating MTT: MTT = CBV/CBF')
[handles.mtt_BzD_cbv, handles.mtt_oSVD_cbv] = mtt_with_CBV_gui(handles); %seconds

%MTT as area under residue function________________________________________
waitbar(3/7, handles.p_bar, 'Calculating MTT: {MTT = \int R(t)}')
[handles.mtt_BzD_no_cbv, handles.mtt_oSVD_no_cbv] = mtt_no_CBV_gui(handles);%seconds

%CBV as CBF*MTT____________________________________________________________
waitbar(4/7, handles.p_bar, 'Calculating CBV: {CBV = CBF \cdot MTT}')
[handles.cbv_BzD_mtt, handles.cbv_oSVD_mtt] = cbv_using_mtt_gui(handles); %ml/100g

%TTP from input data, and from fitted data_________________________________
if (1)
    [handles.ttp_from_data, handles.ttp_fit_BzD, handles.ttp_fit_oSVD] = find_ttp_gui(handles); %seconds
end

%R10 and R50_______________________________________________________________
if handles.calculate_r10 || handles.calculate_r50
    waitbar(5/7, handles.p_bar, 'Calculating R10 and/or R50')
    [handles.r10_BzD, handles.r50_BzD, handles.r10_svd, handles.r50_svd] = r10_r50_gui(handles); %in seconds
end

%OEF and CMRO2_____________________________________________________________
guidata(hObject, handles)
if handles.calculate_oef || handles.calculate_cmro2
    waitbar(6/7, handles.p_bar, 'Calculating OEF and/or CMRO2: User input required')
    [handles.oef_BzD, handles.cmro2_BzD, handles.oef_svd, handles.cmro2_svd] = oef_cmro2_gui(handles);
end
set(handles.slider1, 'Callback', {@slider1_Callback, hObject})
set(handles.slider8, 'Callback', {@slider8_Callback, hObject})

%Done______________________________________________________________________
waitbar(1, handles.p_bar, 'Done.')
delete(handles.p_bar)
set(handles.edit7, 'String', 'Analysis is complete.' ,'ForeGroundColor', 'blue')

%% DISPLAY RESULTS
    handles = ready_the_axes(handles);
    guidata(hObject, handles)
    if handles.show_cbf
        if handles.BzD; handles.cbf_data =  handles.fitd_cbf; else; handles.cbf_data =  handles.fitd_cbf_svd; end
        set(handles.cbf_image, 'CData', imrotate(handles.cbf_data(:,:,1), 90))
        title(handles.cbf_axes, 'CBF [ml/100g/min] : Slice number 1')
        set(handles.cbf_slider, 'Visible', 'on')
        set(handles.cbf_slider, 'Min', 1, 'Max', handles.img_size(3), ...
            'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', 1)
        clim_max = double(max(handles.cbf_data(isfinite(handles.cbf_data)))) + 1E-10;
        if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/999; end
        if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end
        set(handles.cbf_window_slider_up, 'Visible', 'on') %color setting
        set(handles.cbf_window_slider_up, 'Min', 10E-10, 'Max', clim_max, ...
            'SliderStep', [step_small, step_big], 'Value', clim_max)
        set(handles.cbf_window_slider_low, 'Visible', 'on') %color setting
        set(handles.cbf_window_slider_low, 'Min', 0, 'Max', clim_max-1E-10, ...
            'SliderStep', [step_small, step_big], 'Value', 0)
        guidata(hObject, handles)
    end
    if handles.show_mtt
        if handles.BzD; handles.mtt_data =  handles.mtt_BzD_no_cbv; else; handles.mtt_data =  handles.mtt_oSVD_no_cbv; end
        set(handles.mtt_image, 'CData', imrotate(handles.mtt_data(:,:,1), 90))
        title(handles.mtt_axes, 'MTT [s] : Slice number 1')
        set(handles.mtt_slider, 'Visible', 'on')
        set(handles.mtt_slider, 'Min', 1, 'Max', handles.img_size(3), ...
            'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', 1)
        clim_max = double(max(handles.mtt_data(isfinite(handles.mtt_data)))) + 1E-10;
        if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/999; end
        if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end
        set(handles.mtt_window_slider_up, 'Visible', 'on') %color setting
        set(handles.mtt_window_slider_up, 'Min', 10E-10, 'Max', clim_max, ...
            'SliderStep', [step_small, step_big], 'Value', clim_max)
        set(handles.mtt_window_slider_low, 'Visible', 'on') %color setting
        set(handles.mtt_window_slider_low, 'Min', 0, 'Max', clim_max-1E-10, ...
            'SliderStep', [step_small, step_big], 'Value', 0)
        guidata(hObject, handles)
    end


%% CHECK CORRELATIONS (REDUNDANT)
% corerlation between delay and ttp
if (0)
    make_plots = true;
    slice = 9;
    [rho_ttp_delay_BzD, p_ttp_delay_BzD] = correlate_gui(ttp_fit_BzD, fitd_delay, mask, slice, make_plots);
    [rho_ttp_delay_oSVD, p_ttp_delay_oSVD] = correlate_gui(ttp_fit_oSVD, fitd_delay_svd, mask, slice, make_plots);
    
    %correlation between mtt and delay
    [rho_mtt_delay_BzD, p_mtt_delay_BzD] = correlate_gui(mtt_BzD_no_cbv, fitd_delay, mask, slice, make_plots);
    [rho_mtt_delay_oSVD, p_mtt_delay_oSVD] = correlate_gui(mtt_oSVD_no_cbv, fitd_delay_svd, mask, slice, make_plots);
end

childlist = get(handles.figure1, 'Children');
for c = 1:numel(childlist)
    try set(childlist(c), 'Enable', 'on'); catch; end
end
% set(handles.figure1, 'pointer', 'arrow')
guidata(hObject, handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)


function edit7_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement ON PLOT AXES 1. Changes slice
function slider1_Callback(hObject, ~, ~)
handles = guidata(hObject);
index = round(get(hObject, 'Value'));       
                if handles.show_cbf && handles.cbf_slider == handles.slider1
                    set(handles.cbf_image, 'CData', imrotate(handles.cbf_data(:,:,index), 90) );
                    title(handles.axes1, ['CBF [ml/100g/min] : Slice number ' num2str(index)])
                end
                if handles.show_mtt && handles.mtt_slider == handles.slider1
                    set(handles.mtt_image, 'CData', imrotate(handles.mtt_data(:,:,index), 90) );
                    title(handles.axes1, ['MTT [s] : Slice number ' num2str(index)])
                end

guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement ON PLOT AXES 2. Changes slice
function slider2_Callback(hObject, eventdata, handles)
index = round(get(hObject, 'Value'));
            if handles.show_cbf && handles.cbf_slider == handles.slider2
                set(handles.cbf_image, 'CData', imrotate(handles.cbf_data(:,:,index), 90) );
                title(handles.axes2, ['CBF [ml/100g/min] : Slice number ' num2str(index)])
            end
            if handles.show_mtt && handles.mtt_slider == handles.slider2
                set(handles.mtt_image, 'CData', imrotate(handles.mtt_data(:,:,index), 90) );
                title(handles.axes2, ['MTT [s] : Slice number ' num2str(index)])
            end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in DECONVOLUTION ALGORITHM SELECTION.
function popupmenu7_Callback(hObject, eventdata, handles)
algos = cellstr(get(hObject,'String'));
algo = algos{get(hObject,'Value')};
handles.BzD = false; handles.with_delay = false; handles.with_dispersion = false;
handles.do_SVD = false;
handles.do_oSVD = false; handles.do_cSVD = false; handles.do_sSVD = false;
switch algo
    case 'BzD (+ delay)'
        handles.BzD = true;
        handles.with_delay = true;
    case 'BzD'
        handles.BzD = true;
    case 'BzD (+ VTF)'
        handles.BzD = true;
        handles.with_dispersion = true;
    case 'BzD (+ delay + VTF)'
        handles.BzD = true;
        handles.with_delay = true;
        handles.with_dispersion = true;
    case 'oSVD'
        handles.do_oSVD = true; handles.do_SVD = true;
    case 'cSVD'
        handles.do_cSVD = true; handles.do_SVD = true;
    case 'sSVD'
        handles.do_sSVD = true; handles.do_SVD = true;
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% --- Executes on button press in LOAD DATA ONTO AXES1.
function pushbutton8_Callback(hObject, eventdata, handles)
axes1_visualise_data(handles)
%--------------------------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------------------------
% --- Executes on button press in SAVE DATA FROM AXES1.
function pushbutton9_Callback(hObject, eventdata, ~)
    handles = guidata(hObject);
    handles.ext_plot = handles.ext_plot + 1; %allow opening of a limitless number of figures
    handles.external.(['fig' (num2str(handles.ext_plot))]) = figure();
    handles.external.(['ax' (num2str(handles.ext_plot))]) = axes('Parent', handles.external.(['fig' (num2str(handles.ext_plot))]));
    copyobj(get(handles.axes1, 'children'), handles.external.(['ax' (num2str(handles.ext_plot))]))
    colormap( handles.external.(['ax' (num2str(handles.ext_plot))]), colormap(handles.axes1))
    handles.external.(['ax' (num2str(handles.ext_plot))]).CLim = handles.axes1.CLim;
    handles.external.(['ax' (num2str(handles.ext_plot))]).YDir = handles.axes1.YDir;
    handles.external.(['ax' (num2str(handles.ext_plot))]).XLim = handles.axes1.XLim;
    handles.external.(['ax' (num2str(handles.ext_plot))]).YLim = handles.axes1.YLim;
    handles.external.(['ax' (num2str(handles.ext_plot))]).XLabel = handles.axes1.XLabel;
    handles.external.(['ax' (num2str(handles.ext_plot))]).YLabel = handles.axes1.YLabel;
    set(handles.external.(['ax' (num2str(handles.ext_plot))]),  {'XGrid', 'YGrid', 'ZGrid'}, get(handles.axes1, {'XGrid', 'YGrid', 'ZGrid'}))
    title(handles.external.(['ax' (num2str(handles.ext_plot))]), handles.axes1.Title.String, 'Interpreter', 'None')
guidata(hObject, handles)
%------------------------------------------------------------------------------------------------------------------

% --- Executes on button press in LOAD DATA ONTO AXES2.
function pushbutton10_Callback(hObject, eventdata, handles)
axes2_visualise_data(handles)
%---------------------------------------------------------------------------------------------------------

% --- Executes on button press in SAVE DATA FROM AXES2.
function pushbutton11_Callback(hObject, eventdata, handles)
    handles.ext_plot = handles.ext_plot + 1; %allow opening of a limitless number of figures
    handles.external.(['fig' (num2str(handles.ext_plot))]) = figure();
    handles.external.(['ax' (num2str(handles.ext_plot))]) = axes('Parent', handles.external.(['fig' (num2str(handles.ext_plot))]));
    copyobj(get(handles.axes2, 'children'), handles.external.(['ax' (num2str(handles.ext_plot))]))
    colormap( handles.external.(['ax' (num2str(handles.ext_plot))]), colormap(handles.axes2))
    handles.external.(['ax' (num2str(handles.ext_plot))]).CLim = handles.axes2.CLim;
    handles.external.(['ax' (num2str(handles.ext_plot))]).YDir = handles.axes2.YDir;
    handles.external.(['ax' (num2str(handles.ext_plot))]).XLim = handles.axes2.XLim;
    handles.external.(['ax' (num2str(handles.ext_plot))]).YLim = handles.axes2.YLim;
    handles.external.(['ax' (num2str(handles.ext_plot))]).XLabel = handles.axes2.XLabel;
    handles.external.(['ax' (num2str(handles.ext_plot))]).YLabel = handles.axes2.YLabel;
    title(handles.external.(['ax' (num2str(handles.ext_plot))]), handles.axes2.Title.String, 'Interpreter', 'None')
guidata(hObject, handles)

% --- Executes on button press in ADDITIONAL CALCULATIONS: OEF.
function checkbox26_Callback(hObject, eventdata, handles)
handles.calculate_oef = get(hObject,'Value');
handles.calculate_cmro2 = get(hObject,'Value');
guidata(hObject, handles)

% --- Executes on button press in ADDITIONAL CALCULATIONS: R10.
function checkbox28_Callback(hObject, eventdata, handles)
handles.calculate_r10 = get(hObject,'Value');
handles.calculate_r50 = get(hObject,'Value');
guidata(hObject, handles)



% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement CHANGING WINDOW SETTING ON AXES1
function slider6_Callback(hObject, eventdata, handles)
color_max = get(hObject,'Value');
if color_max == 0; color_max = color_max + 1E-10; end
color_min = get(handles.slider10, 'Value');
if color_min > color_max; color_min = color_max-1E-10; set(handles.slider10, 'Value', color_min);end
handles.axes1.CLim = [color_min color_max];
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement CHANGING WINDOW SETTING ON AXES2
function slider7_Callback(hObject, eventdata, handles)
color_max = get(hObject,'Value');
if color_max == 0; color_max = color_max + 1E-10; end
color_min = get(handles.slider11, 'Value');
if color_min > color_max; color_min = color_max-1E-10; set(handles.slider11, 'Value', color_min);end
handles.axes2.CLim = [color_min color_max];
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function slider7_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function uipanel4_CreateFcn(hObject, eventdata, handles)


% --- Executes on slider movement HORIZONTAL SCROLLING ON AXES1.
function slider8_Callback(hObject, ~, ~)
handles= guidata(hObject);
time_point = round(get(hObject,'Value'));
current_slice = round(get(handles.slider1,'Value'));
if handles.plotting_4d_data
    xlim(handles.axes1, [0 time_point])
elseif handles.data_is_for_aif
%     set(handles.aif_conc_image, 'CData', squeeze(handles.aif_conc_data(:,:,current_slice, time_point))  )
%     title(handles.axes1, ['Concentration' ': Slice ' num2str(current_slice) ' : Time-point ' num2str(time_point)],...
%         'Interpreter', 'None')
elseif handles.choosing_region_to_analyse
    slice_indx = round(get(handles.slider1,'Value'));
    set(handles.dsc_data_image, 'CData', squeeze(handles.dsc_data_filtered(:,:,slice_indx,time_point))  )
    title(handles.axes1, ['DSC signal: Slice ' num2str(slice_indx), ' : Time point ', num2str(time_point)])
else
    set(handles.axes1_image, 'CData', squeeze(handles.axes1_data(:,:,current_slice, time_point))  )
    title(handles.axes1, [handles.axes1_filename ': Slice ' num2str(current_slice) ' : Time-point ' num2str(time_point)],...
        'Interpreter', 'None')
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function slider8_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement HORIZONTAL SCROLLING ON AXES2.
function slider9_Callback(hObject, eventdata, handles)
time_point = round(get(hObject,'Value'));
current_slice = round(get(handles.slider2,'Value'));
if handles.plotting_4d_data
    xlim(handles.axes2, [0 time_point])
else
    set(handles.axes2_image, 'CData', squeeze(handles.axes2_data(:,:,current_slice, time_point))  )
    title(handles.axes2, [handles.axes2_filename ': Slice ' num2str(current_slice) ' : Time-point ' num2str(time_point)],...
        'Interpreter', 'None')
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function slider9_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function text27_CreateFcn(hObject, eventdata, handles)


% --- Executes on button press in CHANGE COLORMAP ON AXES1.
function pushbutton15_Callback(hObject, eventdata, handles)
if handles.cmap_index == 6; handles.cmap_index = 0; end
handles.cmap_index = handles.cmap_index + 1;
colormap(handles.axes1, (handles.colormaps{handles.cmap_index}))
drawnow
guidata(hObject, handles)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
if handles.cmap_index == 6; handles.cmap_index = 0; end
handles.cmap_index = handles.cmap_index + 1;
colormap(handles.axes2, (handles.colormaps{handles.cmap_index}))
guidata(hObject, handles)


% --- Executes on button press in COLORBAR ON AXES1.
function togglebutton2_Callback(hObject, eventdata, handles)
if get(hObject,'Value'); colorbar(handles.axes1); else; delete(handles.axes1.Colorbar); end
guidata(hObject, handles)


% --- Executes on button press in togglebutton3.
function togglebutton3_Callback(hObject, eventdata, handles)
if get(hObject,'Value'); colorbar(handles.axes2); else; delete(handles.axes2.Colorbar); end
guidata(hObject, handles)


% --- Executes on slider movement. Lower Window setting on AXES1
function slider10_Callback(hObject, eventdata, handles)
color_max = get(handles.slider6,'Value');
if color_max == 0; color_max = color_max + 1E-10; end
color_min = get(hObject,'Value');
if color_min > color_max; color_min = color_max-1E-10; set(hObject, 'Value', color_min);end
handles.axes1.CLim = [color_min color_max];
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function slider10_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement. Lower Window setting on AXES2
function slider11_Callback(hObject, eventdata, handles)
color_max = get(handles.slider7,'Value');
if color_max == 0; color_max = color_max + 1E-10; end
color_min = get(hObject,'Value');
if color_min > color_max; color_min = color_max-1E-10; set(hObject, 'Value', color_min);end
handles.axes2.CLim = [color_min color_max];
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function slider11_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on key press with focus on pushbutton4 and none of its controls.
function pushbutton4_KeyPressFcn(hObject, eventdata, handles)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over pushbutton4.
function pushbutton4_ButtonDownFcn(hObject, eventdata, handles)



% --- Executes on button press in pushbutton17. ROTATE IMAGE ON AXES1
function pushbutton17_Callback(hObject, eventdata, handles)
try
    axes1_ob = findobj(handles.axes1, 'Type', 'image');
    axes1_ob.CData = imrotate(axes1_ob.CData, -90);
catch
end
guidata(hObject, handles)

% --- Executes on button press in pushbutton18. ROTATE IMAGE ON AXES1
function pushbutton18_Callback(hObject, eventdata, handles)
try
    axes2_ob = findobj(handles.axes2, 'Type', 'image');
    axes2_ob.CData = imrotate(axes2_ob.CData, -90);
catch
end
guidata(hObject, handles)


% --- Executes on button press in pushbutton19.

% hObject    handle function pushbutton19_Callback(hObject, eventdata, handles)to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function pushbutton19_Callback(hObject, eventdata, handles)
if isfield(handles, 'file_folder')
[file, path] = uigetfile( '*.*','Select DSC-MRI data to analyse', handles.file_folder);
else
  [file, path] = uigetfile('*.*','Select DSC-MRI data to analyse');  
end
if file == 0 %if user cancels file selection dialog
    wrap_text = 'No data selected for analysis.'; %report this on display panel
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
    return %terminate execution
end
reset_axes(handles)
path_to_data = fullfile(path, file); %otherwise generate full file path
[file_folder, file_name, file_ext] = fileparts(path_to_data); %extract folder name, file name and file extension
%TODO
%Check that file_ext is nii
%Display file path in a separate panel on main window. This panel
%is static
handles.file_folder = file_folder;
guidata(hObject, handles)
%__________________________________________________________________________
[dsc_data, mask, header_info] = load_dsc_data_gui(path_to_data, handles); %this function loads the data, generates a mask, and returns header info
handles.file_name = file_name;
handles.header_info = header_info;
handles.dsc_data = dsc_data;
handles.mask = mask;
handles.img_size = size(dsc_data);
handles.path_to_data = path_to_data;

if isequal(dsc_data, false); return; end
%display some message on the static panel
%GET TE TR
[te, tr] = get_te_tr(handles); %get te and tr
if te == false && tr == false %if that fails, report and terminate
    wrap_text = 'Analysis terminated by user.';  set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r'); return
else
    handles.te = te; handles.tr = tr;
end
%TODO

handles.dsc_data_loaded = true; 
wrap_text = {'- - - - DSC data path - - - -'; handles.path_to_data}; %report this on display panel
set(handles.edit14, 'String', wrap_text, 'ForegroundColor', 'b')

%Target folder for saving results
% handles.target_folder = generate_target_folder(handles);

% GENERATE WHOLE-BRAIN SIGNAL CURVE
handles.baseline_index = 4:10; %baseline images: first 10 images excluding the first 3
handles.tail_index = handles.img_size(4)-5:handles.img_size(4); %tail: last five images
[~, t_min_signal] = whole_brain_curve_gui(handles); %generate it
handles.t_min_signal = t_min_signal; %save the time point where minimum signal occurs
report('DSC-MRI data loaded successfully.', handles)
guidata(hObject, handles)




% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% GET AIF
    if ~handles.dsc_data_loaded
        wrap_text = 'AIF selection requires DSC data. Load data to proceed.'; %check that data exists
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r') %inform user
        return
    end
    
    %undo restrict analysis
    if isfield(handles, 'old_mask'); handles.mask = handles.old_mask; end
    
set([handles.pushbutton31, handles.pushbutton32,handles.pushbutton33, handles.pushbutton34, handles.pushbutton35], 'Visible', 'off') %initially turn off all these buttons
set([handles.pushbutton24, handles.pushbutton25,handles.pushbutton26, handles.pushbutton27, handles.pushbutton30], 'Visible', 'off') %initially turn off all these buttons
set(handles.figure1, 'WindowButtonMotionFcn', []);
reset_axes(handles, handles.axes2)
if handles.calculate_aif %if user wants to select AIF from input data 
    %remove these buttons as they may cause confusion
    set([handles.pushbutton8, handles.pushbutton9, handles.pushbutton10, handles.pushbutton11], 'Visible', 'off')
    handles.roi_counter = 0; %needed to keep track of how many rois have been drawn 
    %reset the axes
    cla(handles.axes1, 'reset')
    cla(handles.axes2, 'reset')
    set(handles.axes1,'XTick', [], 'YTick', [])
    set(handles.axes2,'XTick', [], 'YTick', [])
    set(handles.axes1, 'box', 'on')
    set(handles.axes2, 'box', 'on')
    %______________________________________________________________________
    %turn on all image-related  objects on axes1; this is where the
    %concentration image will be shown for the user to place a ROI. Turn
    %off all image-related objects on axes2; this is where the AIFs will be
    %plotted
    set(handles.slider1,'Visible', 'on');
    set(handles.slider2,'Visible', 'off');
    set(handles.slider6, 'Visible', 'on');
    set(handles.slider10, 'Visible', 'on');
    set(handles.pushbutton15, 'Visible', 'on')
    set(handles.pushbutton16, 'Visible', 'off')
    set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
    set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes2
    set(handles.togglebutton2, 'Visible', 'on')
    set(handles.togglebutton3, 'Visible', 'off')
    %______________________________________________________________________
    wrap_text = 'Preparing for AIF selection ... '; %it will take a short while to prepare the concentration image so
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b') %inform user
    opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
    handles.wait_for_data = msgbox('\fontsize{10}Preparing for AIF selection. Please wait...', 'AIF selection','none' , opts);
    aif_handles.aif_slice = ceil(handles.img_size(3)/2); %choose middle slice as initial slice for AIF selection; can be changed using slider
    aif_handles.aif_time_point = handles.t_min_signal; %save time point for minimum signal
    handles.aif_handles = aif_handles; %aif_handles contains settings for AIF selection
    guidata(hObject, handles) %update app data
    handles.hObject = hObject; %supply AIF selection algorithm with handle to current object (pushbutton4)
    if handles.calculate_aif_auto %if user chose the automatic AIF selector
        find_aif_auto_gui(handles); %call corresponding function
    end
    if handles.calculate_aif_semi_auto %user chose semi-automatic
        find_aif_semi_auto_gui(handles); %call corresponding function
    end
    if handles.calculate_aif_manual %user wants to manually select AIF
        find_aif_gui(handles); %call responsible function
    end
    waitfor(handles.pushbutton30, 'Visible')
    handles = guidata(hObject);
    if ~isfield(handles, 'best_aif')
        report('No AIF selected.', handles, 'b')
        return
    end
    mean_aif_c = handles.best_aif;
    aif_area = handles.aif_area;
    handles.aif_handles = []; %otherwise, clear aif_handles
    wrap_text = 'AIF selection complete.'; %display confirmation of completed AIF selection
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
else %user instead wants to load an existing AIF from file
    %turn off all image-related objects on both axes
    set(handles.slider1,'Visible', 'off');
    set(handles.slider2,'Visible', 'off');
    set(handles.slider6, 'Visible', 'off');
    set(handles.slider10, 'Visible', 'off');
    set(handles.pushbutton15, 'Visible', 'off')
    set(handles.pushbutton16, 'Visible', 'off')
    set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
    set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes2
    set(handles.togglebutton2, 'Visible', 'off')
    set(handles.togglebutton3, 'Visible', 'off')
    %______________________________________________________________________
    wrap_text = 'Select AIF';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    [file, path] = uigetfile({'*.*'},'Select AIF for analysis.', handles.file_folder); %invoke file-select dialog with this title
    if file == 0 %if user cancels dialog
        wrap_text = 'No AIF selected.' ; %report and terminate execution
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
        return
    end
    path_to_aif = fullfile(path, file); %otherwise generate full path to AIF file
    mean_aif_c = read_this_file(path_to_aif, 1); %call the function that reads files; output must be 1-dimensional; function raises error otherwise
    if true; t = 0:handles.tr:(handles.img_size(4)-1)*handles.tr; %if you must display aif, generate time-vector
        plot(t, mean_aif_c, 'k-', 'Parent', handles.axes2); %plot
        xlabel( handles.axes2,'t [s]'); ylabel( handles.axes2, '{C_a}(t)'); title( handles.axes2, 'AIF')
    end
end
handles.mean_aif_c = mean_aif_c; %save AIF to handles structure
handles.data_is_for_aif = false; %alert objects that AIF selection is done
handles.aif_exists = true;
tmp = get(handles.edit14, 'String');
tmp{end+1} = '- - - - AIF path - - - -'; 
if handles.calculate_aif
    if isfield(handles, 'path_to_aif'); path_to_aif = handles.path_to_aif; 
    else
       path_to_aif = []; 
    end
end
tmp{end+1} = path_to_aif; 
set(handles.edit14, 'String',tmp)
report('AIF selection successful.', handles)
set(handles.slider1, 'Callback', {@slider1_Callback, hObject})
set(handles.slider8, 'Callback', {@slider8_Callback, hObject})
set([handles.pushbutton8, handles.pushbutton9, handles.pushbutton10, handles.pushbutton11], 'Visible', 'on')
guidata(hObject, handles)




% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% GET VOF
    if ~handles.dsc_data_loaded
        wrap_text = 'VOF selection requires DSC data. Load data to proceed.'; %check that data exists
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r') %inform user
        return
    end
    
        %undo restrict analysis
    if isfield(handles, 'old_mask'); handles.mask = handles.old_mask; end
    
set([handles.pushbutton31, handles.pushbutton32,handles.pushbutton33, handles.pushbutton34, handles.pushbutton35], 'Visible', 'off') %initially turn off all these buttons
set([handles.pushbutton24, handles.pushbutton25,handles.pushbutton26, handles.pushbutton27, handles.pushbutton30], 'Visible', 'off') %initially turn off all these buttons
set(handles.figure1, 'WindowButtonMotionFcn', []);
reset_axes(handles, handles.axes2)
if handles.calculate_vof %if user wants to select VOF from input data 
    set([handles.pushbutton8, handles.pushbutton9, handles.pushbutton10, handles.pushbutton11], 'Visible', 'off')
    handles.roi_counter = 0; %needed to keep track of how many rois have been drawn 
    %reset the axes
    cla(handles.axes1, 'reset')
    cla(handles.axes2, 'reset')
    set(handles.axes1,'XTick', [], 'YTick', [])
    set(handles.axes2,'XTick', [], 'YTick', [])
    set(handles.axes1, 'box', 'on')
    set(handles.axes2, 'box', 'on')
    %______________________________________________________________________
    %turn on all image-related  objects on axes1; this is where the
    %concentration image will be shown for the user to place a ROI. Turn
    %off all image-related objects on axes2; this is where the AIFs will be
    %plotted
    set(handles.slider1,'Visible', 'on');
    set(handles.slider2,'Visible', 'off');
    set(handles.slider6, 'Visible', 'on');
    set(handles.slider10, 'Visible', 'on');
    set(handles.pushbutton15, 'Visible', 'on')
    set(handles.pushbutton16, 'Visible', 'off')
    set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
    set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes2
    set(handles.togglebutton2, 'Visible', 'on')
    set(handles.togglebutton3, 'Visible', 'off')
    %______________________________________________________________________
    wrap_text = 'Preparing for VOF selection ... '; %it will take a short while to prepare the concentration image so
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b') %inform user
    opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
    handles.wait_for_data = msgbox('\fontsize{10}Preparing for VOF selection. Please wait...', 'VOF selection','none' , opts);
    vof_handles.vof_slice = ceil(handles.img_size(3)/2); %choose middle slice as initial slice for VOF selection; can be changed using slider
    vof_handles.vof_time_point = handles.t_min_signal; %save time point for minimum signal
    handles.vof_handles = vof_handles; %aif_handles contains settings for VOF selection
    guidata(hObject, handles) %update app data
    handles.hObject = hObject; %supply VOF selection algorithm with handle to current object (pushbutton4)
    if handles.calculate_vof_auto %if user chose the automatic VOF selector
        find_vof_auto_gui(handles); %call corresponding function
    end
    if handles.calculate_vof_semi_auto %user chose semi-automatic
        find_vof_semi_auto_gui(handles); %call corresponding function
    end
    if handles.calculate_vof_manual %user wants to manually select VOF
        find_vof_gui(handles); %call responsible function
    end
    waitfor(handles.pushbutton30, 'Visible')
    handles = guidata(hObject);
        if ~isfield(handles, 'best_vof')
        report('No VOF selected.', handles, 'b')
        return
    end
    mean_vof_c = handles.best_vof;
    vof_area = handles.vof_area;
    handles.vof_handles = []; %otherwise, clear vof_handles
    wrap_text = 'VOF selection complete.'; %display confirmation of completed VOF selection
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
else %user instead wants to load an existing VOF from file
    %turn off all image-related objects on both axes
    set(handles.slider1,'Visible', 'off');
    set(handles.slider2,'Visible', 'off');
    set(handles.slider6, 'Visible', 'off');
    set(handles.slider10, 'Visible', 'off');
    set(handles.pushbutton15, 'Visible', 'off')
    set(handles.pushbutton16, 'Visible', 'off')
    set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
    set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes2
    set(handles.togglebutton2, 'Visible', 'off')
    set(handles.togglebutton3, 'Visible', 'off')
    %______________________________________________________________________
    wrap_text = 'Select VOF';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    [file, path] = uigetfile({'*.*'},'Select VOF for analysis.', handles.file_folder); %invoke file-select dialog with this title
    if file == 0 %if user cancels dialog
        wrap_text = 'No VOF selected. Analysis terminated.'; %report and terminate execution
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
        return
    end
    path_to_vof = fullfile(path, file); %otherwise generate full path to VOF file
    mean_vof_c = read_this_file(path_to_vof, 1); %call the function that reads files; output must be 1-dimensional; function raises error otherwise
    if true; t = 0:handles.tr:(handles.img_size(4)-1)*handles.tr; %if you must display aif, generate time-vector
        plot(t, mean_vof_c, 'k-', 'Parent', handles.axes2); %plot
        xlabel( handles.axes2,'t [s]'); ylabel( handles.axes2, '{C_a}(t)'); title( handles.axes2, 'VOF')
    end
end
handles.mean_vof_c = mean_vof_c; %save VOF to handles structure
handles.data_is_for_vof = false; %alert objects that VOF selection is done
handles.vof_exists = true;
tmp = get(handles.edit14, 'String');
tmp{end+1} = '- - - - VOF path - - - -'; 
if handles.calculate_vof
    if isfield(handles, 'path_to_vof'); path_to_vof = handles.path_to_vof; 
    else
       path_to_vof = []; 
    end
end
tmp{end+1} = path_to_vof; 
set(handles.edit14, 'String',tmp)
report('VOF selection successful.', handles)
set(handles.slider1, 'Callback', {@slider1_Callback, hObject})
set(handles.slider8, 'Callback', {@slider8_Callback, hObject})
set([handles.pushbutton8, handles.pushbutton9, handles.pushbutton10, handles.pushbutton11], 'Visible', 'on')
guidata(hObject, handles)

% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oef_cmro2_gui_indep(handles)
guidata(hObject, handles)

% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% r10_r50_gui_indep(handles)
report('Independent calculation of R10/R50 is yet to be implemented.', handles, 'r')
guidata(hObject, handles)

% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton26. Add ROI
function pushbutton26_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%This callback will be implemented by functions that need it



% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit13_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit14_Callback(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit14 as text
%        str2double(get(hObject,'String')) returns contents of edit14 as a double


% --- Executes during object creation, after setting all properties.
function edit14_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton27.
function pushbutton27_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton30.
function pushbutton30_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton31.
function pushbutton31_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton32.
function pushbutton32_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton33.
function pushbutton33_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton34.
function pushbutton34_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton35.
function pushbutton35_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
