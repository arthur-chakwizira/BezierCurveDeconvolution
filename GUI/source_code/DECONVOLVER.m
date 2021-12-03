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
%      arthurchakiwizira@gmail.com
%      Medical Radiation Physics, Lund University, Sweden
%
%      Part of the code base is an extension of preliminary work by:
%      Andr? Ahlgren

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
set(handles.pushbutton14, 'Visible', 'off') %Continue button; this is used by one of the AIF/VOF selection functions
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
handles.choosing_region_to_analyse = false; %alert objects that user is placing ROI to restrict analysis to certain region(s)
handles.data_is_for_aif = false; %alert objects that user is currently selecting AIF
handles.data_is_for_vof = false; %alert objects that user is currently selecting VOF
%all the variables below are used to check if some function has reserved a given slider
handles.cbf_slider = 0;
handles.mtt_slider = 0;
handles.delay_slider = 0;
handles.rt_slider = 0;
handles.oef_slider = 0;
handles.cmro2_slider = 0;
handles.r10_slider = 0;
handles.r50_slider = 0;
%_________________________________________________________________________
set(handles.axes1,'XTick', [], 'YTick', []) %set default appearance for axes
set(handles.axes2,'XTick', [], 'YTick', [])
set(handles.axes1, 'box', 'on')
set(handles.axes2, 'box', 'on')
%__________________________________________________________________________
handles.load_existing_data = true; %loading data not simulating data
handles.simulate_data = false; %simulating data; will invoke a DSC simulation program in a future release
set(handles.checkbox24, 'value', 1) %edit_data_setup;  allow user to choose whether to edit mask settings, etc
handles.kappa = 0.705; %kappa; the hematocrit constant
handles.calculate_aif = false; %will invoke AIF selection algorithm if true
set(handles.checkbox7, 'value', 1) %display_aif
handles.include_vof = false; %choose whether to include a VOF for PVE correction
set(handles.checkbox8, 'value', 0) %display_vof
handles.BzD = true; handles.with_delay = true; handles.with_dispersion = false; %initial deconvolution algorithm settings
handles.do_SVD = false; handles.do_oSVD = false; handles.do_cSVD = false; handles.do_sSVD = false; %deconvolution algorithm settings
set(handles.checkbox25, 'value', 1) %edit_algo_setup; let user choose to edit settings such as kappa, Oscillation Index, etc
set(handles.checkbox26, 'value', 1) %calculate_oef
set(handles.checkbox27, 'value', 1) %calculate_cmro2
set(handles.checkbox28, 'value', 0) %calculate_r10
set(handles.checkbox29, 'value', 0) %calculate_r50
handles.post_display = true; %always true; means that results are diisplayed after analysis
handles.live_display = false; %if true, results (some) will be shown during deconvolution
%the following are initial settings for what is displayed on the axes; Default is CBF on left axes and MTT on right axes)
handles.show_cbf = true; handles.cbf_axes = handles.axes1; handles.cbf_slider = handles.slider1; handles.cbf_window_slider_up = handles.slider6; handles.cbf_window_slider_low = handles.slider10;
handles.show_mtt = true; handles.mtt_axes = handles.axes2; handles.mtt_slider = handles.slider2; handles.mtt_window_slider_up = handles.slider7; handles.mtt_window_slider_low = handles.slider11;
handles.show_oef = false;handles.show_cmro2 = false;
handles.show_r10 = false; handles.show_r50 = false;
handles.show_delay = false; handles.plot_residue_funcs = false;
%__________________________________________________________________________
set(handles.checkbox15, 'value', 1) %save_results
handles.save_format = 'NIFTI (*.nii)'; %save_format, default is NIFTI
set(handles.checkbox14, 'value', 1) %notify_when_done, PC will beep if true
wrap_text = 'ANALYSIS OF DSC-MRI DATA.   VERIFY CONFIGURATION AND RUN TO SELECT DATA.'; %default text shown on display panel
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
% Choose default command line output for DECONVOLVER
handles.output = hObject;
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


% --- Executes on button press in DISPLAY AIF CHECKBOX.
function checkbox7_Callback(hObject, eventdata, handles)
handles.display_aif = get(hObject,'Value');
guidata(hObject, handles)


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


% --- Executes on button press in DISPLAY VOF CHECKBOX.
function checkbox8_Callback(hObject, eventdata, handles)
handles.display_vof = get(hObject,'Value');
guidata(hObject, handles)


% --- Executes on button press in NOTIFY WHEN DONE.
function checkbox14_Callback(hObject, eventdata, handles)
handles.notify_when_done = get(hObject,'Value');
guidata(hObject, handles)


% --- Executes on button press in SAVE RESULTS.
function checkbox15_Callback(hObject, eventdata, handles)
handles.save_results = get(hObject,'Value');
guidata(hObject, handles)


% --- Executes when user presses RUN BUTTON.
function pushbutton4_Callback(hObject, eventdata, handles)
%THIS IS THE MAIN PART OF THE PROGRAM

%% ANALYSIS OF DSC-MRI DATA

%% GET CONFIGURATION SETTINGS <>
set(handles.pushbutton14, 'Visible', 'off') %button is used by Manual AIF/VOF selection algorithm
handles.data_is_for_aif = false;%if true, slider1 will assume that concentration data is being shown on axes1
handles.choosing_region_to_analyse = false; %if true, slider1 will assume that dsc_data is being shown on axes1
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
handles.edit_data_setup = get(handles.checkbox24,'Value'); %edit input data setup (mask, mask threshold, etc)
handles.display_aif = get(handles.checkbox7,'Value'); %display_aif
handles.display_vof = get(handles.checkbox8,'Value'); %display_vof
handles.edit_algo_setup = get(handles.checkbox25,'Value'); %edit algorithm setup (hematocrit kappa, OI, Psvd, etc)
handles.notify_when_done = get(handles.checkbox14,'Value'); %notify_when_done (beep when done)
handles.save_results = get(handles.checkbox15,'Value'); %save_results
handles.calculate_oef = get(handles.checkbox26, 'Value'); %calculate_oef
handles.calculate_cmro2 = get(handles.checkbox27, 'Value'); %calculate_cmro2
handles.calculate_r10 = get(handles.checkbox28, 'Value'); %calculate_r10
handles.calculate_r50 = get(handles.checkbox29, 'Value'); %calculate_r50
%restore axes to their state when window opened
cla(handles.axes1, 'reset')
cla(handles.axes2, 'reset')
set(handles.axes1,'XTick', [], 'YTick', [])
set(handles.axes2,'XTick', [], 'YTick', [])
set(handles.axes1, 'box', 'on')
set(handles.axes2, 'box', 'on')
%__________________________________________________________________________
%check if selected file save format is optimal/feasible
if ~strcmp(handles.save_format, 'NIFTI (*.nii)') && ~strcmp(handles.save_format, 'DICOM (*.dcm)') %user has not selected NIFTI or Dicom
    if strcmp(handles.save_format, 'JPEG Image (*.jpg)') || strcmp(handles.save_format, 'Portable Network Graphics (*.png)') %forbid jpg or png save formats
        opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
        errordlg(['\fontsize{10}The selected file save format: (' handles.save_format ') is not appropriate for saving of results of DSC-MRI analysis.'],...
            'Inappropriate file save format', opts); %raise error
        wrap_text = 'Analysis terminated due to inappropriate file save format. Change format and re-run to select data.';
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
        return %terminate execution
    else  %otherwise warn that the selected format is not optimal
        opts.Interpreter = 'tex';
        opts.Default = 'No';
        response = questdlg(['\fontsize{10}The currently selected file save format: (' handles.save_format ') is not recommended for saving of data with' ...
            ' dimensions higher than 2. Best performance is obtained with format NIFTI (*.nii). Proceed with current format?'], ...
            'Save format warning','No', 'Yes', opts);
        switch response
            case 'Yes'
            case 'No'
                wrap_text = 'Analysis terminated by user. Re-run to select data.';
                set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                return
        end
    end
end
%__________________________________________________________________________
%Warn that saving to Dicom is in the current version unreliable and forbid
%execution
if strcmp(handles.save_format, 'DICOM (*.dcm)')
    opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
    errordlg(['\fontsize{10}The selected file save format: (' handles.save_format ') is currently unreliable for saving of results of DSC-MRI analysis.'],...
        'Unreliable file save format', opts); %raise error
    wrap_text = 'Analysis terminated due to unreliable file save format. Change format and re-run to select data.';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
    return %terminate execution
end

%Warn that the quantites OEF, CMRO2, R10 and R50 can only be displayed
%after analysis
if handles.live_display && (handles.show_oef || handles.show_cmro2 || handles.show_r10 || handles.show_r50)
    opts.Interpreter = 'tex';
    opts.Default = 'Yes';
    response = questdlg(['\fontsize{10} The quantities OEF, CMRO2, R10 and R50 will only be displayed when analysis is complete.'...
        ' Proceed?'], ...
        'Results display warning','Yes', 'No', opts);
    switch response
        case 'Yes'
            handles.live_display = false;
            guidata(hObject, handles)
        case {'No', ''}
            wrap_text = 'Analysis terminated by user. Re-run to select data.';
            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
            return
    end
end
%__________________________________________________________________________
%Warn that live plotting of residue functions will slow down the program
if handles.live_display && handles.plot_residue_funcs
    opts.Interpreter = 'tex';
    opts.Default = 'After analysis';
    response = questdlg(['\fontsize{10} You have chosen to display residue functions during data analysis. This will drastically slow down'...
        ' deconvolution and is strongly unrecommended. Select display mode.'], ...
        'Results display warning','After analysis', 'During analysis','Terminate analysis', opts);
    switch response
        case 'After analysis'
            handles.live_display = false;
            set(handles.popupmenu8, 'Value', 1)
            guidata(hObject, handles)
        case 'During analysis'
        case {'', 'Terminate analysis'}
            wrap_text = 'Analysis terminated by user. Re-run to select data.';
            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
            return
    end
end
%__________________________________________________________________________

%Select file filter to be used for saving of files (needed for AIF / VOF)
switch handles.save_format
    case 'NIFTI (*.nii)'; file_filter = {'NIFTI *.nii'; 'DICOM *.dcm' ; 'MATLAB *.mat'; 'Text *.txt'; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat' };
    case 'DICOM (*.dcm)'; file_filter = {'DICOM *.dcm'; 'NIFTI *.nii' ; 'MATLAB *.mat'; 'Text *.txt'; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat'};
    case 'MATLAB (*.mat)'; file_filter = {'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii' ; 'Text *.txt' ; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat'};
    case 'Text (*.txt)'; file_filter = {'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii'; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat' };
    case 'DAT (*.dat)'; file_filter = {'DAT *.dat'; 'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii';  'Comma Separated Values *.csv'; 'Excel *.xlsx' };
    case 'Comma Separated Values (*.csv)'; file_filter = {'Comma Separated Values *.csv'; 'DAT *.dat'; 'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii'; 'Excel *.xlsx' };
    case 'Excel (*.xlsx)';  file_filter = {'Excel *.xlsx'; 'Comma Separated Values *.csv'; 'DAT *.dat'; 'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii' };
end
handles.file_filter = file_filter;
%% SIMULATE DSC DATA
%If user chooses to select data, notify that this functionality is not yet
%available
if handles.simulate_data
    opts.Interpreter = 'tex'; opts.WindowStyle = 'non-modal';
    wrap_text = ['\fontsize{10} Data simulation has not yet been implemented in this program. Please load existing DSC-MRI data. For help, contact: ', newline, ...
        'Arthur Chakwizira ',newline, 'arthurchakwizira@gmail.com ',newline, 'Lund University, Sweden'];
    msgbox(wrap_text, 'DSC-MRI data simulation' ,'help', opts);
    guidata(hObject, handles)
    wrap_text = 'Analysis terminated following attempt to call non-existent function. Re-run to select data.';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
    return
end

%% LOAD EXISTING DSC DATA
if handles.load_existing_data %if user chose to load existing data
    wants_new_data = true;
    if isfile('existing_data.mat') %check if there is already data saved
        old_data = load('existing_data.mat');
        opts.Interpreter = 'tex';
        opts.Default = 'Yes';
        prompt = '\fontsize{10}Data previously loaded from the path displayed on the main window already exists. Proceed with this data?';
        wrap_text = old_data.path_to_data; %report this on display panel
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
        response = questdlg(prompt, ... %inform that program will ask for cell range
            'Loading DSC-MRI data','Yes', 'No', 'Cancel', opts);
        switch response
            case 'No'; wants_new_data = true; %user wants to load new data
            case 'Yes'; wants_new_data = false;
            case {'Cancel', ''};     wrap_text = 'Analysis terminated by user. Re-run to select data.'; %invoke file selection dialog with this title
                set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r') %display same message on display panelreturn
                return
        end
    else
        wants_new_data = true;
    end
    
    
    if wants_new_data
        [file, path] = uigetfile('*.*','Select DSC-MRI data to analyse');
        if file == 0 %if user cancels file selection dialog
            wrap_text = 'No data selected for analysis. Re-run to select data.'; %report this on display panel
            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
            return %terminate execution
        end
        path_to_data = fullfile(path, file); %otherwise generate full file path
        [~, ~, file_ext] = fileparts(path_to_data); %extract folder name, file name and file extension
        %Forbid the reading of dicom files as input for DSC analysis
        %This functionality needs more work to be reliable
        %__________________________________________________________________________
        if strcmp(file_ext, '.dcm')
            opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
            dicom_msg = msgbox(['\fontsize{10}Construction of a 4D volume from a series of dicom images may take a while in the current version.'...
                ' Use of the NIFTI format is strongly recommended.'], 'Dicom file format','warn' , opts);
            %         return
            uiwait
            if ~isgraphics(dicom_msg); uiresume; end
        end
        %__________________________________________________________________________
        %The code below displays an input dialog and promts the user to choose to
        %generate a mask, with what threshold, generate whole-brain-curve, display
        %whole-brain-curve
        handles.path_to_data = path_to_data; %save path to loaded data in handles structure
        wrap_text = 'Waiting for user input...';
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
        opts.Interpreter = 'tex';
        prompt = { '\color{blue} \fontsize{10}Make mask (yes/ no)',...
            '\color{blue} \fontsize{10}Mask threshold:', '\color{blue} \fontsize{10} Save mask (yes/ no)', ...
            '\color{blue} \fontsize{10}Create whole-brain signal curve: (yes/ no)',...
            '\color{blue} \fontsize{10}Display whole-brain signal curve: (yes/ no)'};
        dlgtitle = 'DECONVOLVER: Edit input data setup';
        dims = [1 80];
        try  %check if user has entered these choices before
            current_make_mask = num2str(handles.make_mask); if current_make_mask; current_make_mask = 'yes'; else; current_make_mask = 'no'; end
            current_mask_threshold = num2str(handles.mask_threshold); current_save_mask = num2str(handles.save_mask);
            current_create_whole_brain_curve = num2str(handles.create_whole_brain_curve); if current_create_whole_brain_curve; current_create_whole_brain_curve = 'yes'; else; current_create_whole_brain_curve = 'no'; end
            current_make_plots = num2str(handles.make_plots); if current_make_plots; current_make_plots = 'yes'; else; current_make_plots = 'no'; end
        catch %if not, use these default values
            current_make_mask = 'yes';
            current_mask_threshold = '100'; current_save_mask = 'yes';
            current_create_whole_brain_curve = 'yes';
            current_make_plots = 'yes';
        end
        
        defaults = {current_make_mask, current_mask_threshold, current_save_mask,...
            current_create_whole_brain_curve, current_make_plots};
        
        response = inputdlg(prompt, dlgtitle, dims, defaults, opts); %input dialog
        
        if isempty(response) %if user clicks cancel, report and terminate
            wrap_text = 'Analysis terminated by user. Re-run to select data.';
            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
            return
        end
        if strcmp(response{1}, 'no'); handles.make_mask = false; else; handles.make_mask = true; end %extract user inputs
        handles.mask_threshold = str2double(response{2});
        if strcmp(response{3}, 'no'); handles.save_mask = false; else; handles.save_mask = true; end
        if strcmp(response{4}, 'no'); handles.create_whole_brain_curve = false; else; handles.create_whole_brain_curve = true; end
        if strcmp(response{5}, 'no'); handles.make_plots = false; else; handles.make_plots = true; end
        guidata(hObject, handles)
        %__________________________________________________________________________
        [dsc_data, mask, header_info] = load_dsc_data_gui(path_to_data, handles); %this function loads the data, generates a mask, and returns header info
        if isequal(dsc_data, false) && isequal(mask, false); return; end %terminate execution if data loading failed
    else
        header_info = old_data.header_info;
        dsc_data = old_data.dsc_data;
        mask = old_data.mask;
        path_to_data = old_data.path_to_data;
        handles.create_whole_brain_curve = true;
        handles.make_plots = true;
        handles.make_mask = true;
        handles.save_mask = true;
    end
    handles.header_info = header_info;
    handles.dsc_data = dsc_data;
    handles.mask = mask;
    handles.img_size = size(dsc_data);
    handles.path_to_data = path_to_data;
    %__________________________________________________________________________
    %Save loaded data and path to allow quick retrieval on re-run
    save('existing_data.mat', 'dsc_data', 'mask', 'header_info', 'path_to_data');
    %__________________________________________________________________________
    if handles.edit_data_setup
        [te, tr] = get_te_tr(handles); %get te and tr
        if te == false && tr == false %if that fails, report and terminate
            wrap_text = 'Analysis terminated by user. Re-run to select data.';  set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r'); return
        else
            handles.te = te; handles.tr = tr;
        end
    end
    guidata(hObject, handles)
end

%% Target folder for saving results
[file_path, file_name, file_ext] = fileparts(path_to_data);%folder name, file name, file extension
if strcmp(file_ext, '.dcm') %if save format is dicom
    [file_path, file_name, ~] = fileparts(file_path); %get path to folder and folder name
    handles.target_folder = strcat(file_path, '\Results_for_' , file_name); %save results in folder with same name as input folder
else %for other formats
    handles.target_folder = strcat(file_path, '\Results_for_' , file_name); %save in folder with same name asa input file
end
if ~exist(handles.target_folder, 'dir'); mkdir(handles.target_folder); end %if target folder is non-existent, make it

%% GENERATE WHOLE-BRAIN SIGNAL CURVE
handles.baseline_index = 4:10; %baseline images: first 10 images excluding the first 3
handles.tail_index = handles.img_size(4)-5:handles.img_size(4); %tail: last five images
if handles.create_whole_brain_curve %if user wants whole brain dsc curve
    [~, t_min_signal] = whole_brain_curve_gui(handles); %generate it
    handles.t_min_signal = t_min_signal; %save the time point where minimum signal occurs
end
guidata(hObject, handles)

%% GET AIF
if handles.calculate_aif %if user wants to select AIF from input data
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
    wrap_text = 'Loading data for AIF selection ... '; %it will take a short while to prepare the concentration image so
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b') %inform user
    opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
    handles.wait_for_data = msgbox('\fontsize{10}Loading data for AIF selection. Please wait...', 'AIF selection','none' , opts);
    aif_handles.aif_slice = ceil(handles.img_size(3)/2); %choose middle slice as initial slice for AIF selection; can be changed using slider
    aif_handles.display_aif = handles.display_aif; %display AIF, after completed selection
    aif_handles.save_aif = true;
    if ~exist('t_min_signal','var') %if the time point for minimum signal is non-existent (because user chose not to generate whole brain dsc curve)
        current_state = handles.make_plots;
        handles.make_plots = false;
        [~, t_min_signal] = whole_brain_curve_gui(handles); %generate whole brain curve without displaying it
        handles.make_plots = current_state;
    end
    aif_handles.aif_time_point = t_min_signal; %save time point for minimum signal
    handles.aif_handles = aif_handles; %aif_handles contains settings for AIF selection
    guidata(hObject, handles) %update app data
    handles.hObject = hObject; %supply AIF selection algorithm with handle to current object (pushbutton4)
    if handles.calculate_aif_auto %if user chose the automatic AIF selector
        [mean_aif_c, aif_area] = find_aif_auto_gui(handles); %call corresponding function
    end
    if handles.calculate_aif_semi_auto %user chose semi-automatic
        [mean_aif_c, aif_area] = find_aif_semi_auto_gui(handles); %call corresponding function
    end
    if handles.calculate_aif_manual %user wants to manually select AIF
        [mean_aif_c, aif_area] = find_aif_gui(handles); %call responsible function
    end
    if mean_aif_c == false & aif_area == false %if this happens, AIF selection was terminated by user
        set(handles.pushbutton14, 'Visible', 'off')
        wrap_text = 'AIF selection terminated by user.  Rerun to select data. '; %report and terminate execution
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
        return
    end
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
    [file, path] = uigetfile({'*.*'},'Select AIF for analysis.'); %invoke file-select dialog with this title
    if file == 0 %if user cancels dialog
        wrap_text = 'No AIF selected. Analysis terminated.'; %report and terminate execution
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
        return
    end
    path_to_aif = fullfile(path, file); %otherwise generate full path to AIF file
    mean_aif_c = read_this_file(path_to_aif, 1); %call the function that reads files; output must be 1-dimensional; function raises error otherwise
    if handles.display_aif; t = 0:handles.tr:(handles.img_size(4)-1)*handles.tr; %if you must display aif, generate time-vector
        if handles.make_plots; aif_axes = handles.axes2; else; aif_axes = handles.axes1; end %determine which axes to plot on
        plot(t, mean_aif_c, 'k-', 'Parent', aif_axes); %plot
        xlabel(aif_axes,'t [s]'); ylabel(aif_axes, '{C_a}(t)'); title(aif_axes, 'AIF')
    end
end
handles.mean_aif_c = mean_aif_c; %save AIF to handles structure
if ~exist('aif_area', 'var'); t = 0:handles.tr:(handles.img_size(4)-1)*handles.tr; %if area under AIF has not been calculated; generate time-vector
    aif_area = trapz(t, mean_aif_c); end %and calculate area under AIF
handles.aif_area = aif_area; %save area under AIF
handles.data_is_for_aif = false; %alert objects that AIF selection is done

%% GET VOF
if handles.include_vof %if user wants to include VOF (for PVE correction)
    if handles.calculate_vof %and wants to select the VOF from input data
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
        %off all image-related objects on axes2; this is where the VOFs will be
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
        wrap_text = 'Loading data for VOF selection ... '; %it may take a short while to generate concentration images
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b') %inform user
        opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
        handles.wait_for_data = msgbox('\fontsize{10}Loading data for VOF selection. Please wait...', 'VOF selection','none' , opts);
        vof_handles.vof_slice = ceil(handles.img_size(3)/2); %choose middle slice as initial slice for VOF selection; can be changed using slider
        vof_handles.display_vof = true;
        vof_handles.save_vof = true;
        if ~exist('t_min_signal','var') %if the time point for minimum signal does not exist
            current_state = handles.make_plots;
            handles.make_plots = false;
            [~, t_min_signal] = whole_brain_curve_gui(handles); %generate whole brain curve without displaying it
            handles.make_plots = current_state;
        end
        vof_handles.vof_time_point = t_min_signal; %extract time point for minimum signal
        handles.vof_handles = vof_handles; %vof_handles contains settings for VOF selection
        guidata(hObject, handles)
        handles.hObject = hObject;
        if handles.calculate_vof_auto %user chose automatic VOF selection
            [mean_vof_c, vof_area] = find_vof_auto_gui(handles); %call responsible function
        end
        if handles.calculate_vof_semi_auto %user chose semi-automatic VOF selector
            [mean_vof_c, vof_area] = find_vof_semi_auto_gui(handles); %call selector
        end
        if handles.calculate_vof_manual %user chose manual VOF selection
            [mean_vof_c, vof_area] = find_vof_gui(handles); %call responsible function
        end
        if mean_vof_c == false & vof_area == false %if this happens, VOF selection was cancelled by the user
            set(handles.pushbutton14, 'Visible', 'off')
            wrap_text = 'VOF selection terminated by user.  Rerun to select data. '; %report
            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')    %terminate execution
            return
        end
        handles.vof_handles = []; %clear VOF selection handles
        handles.mean_vof_c = mean_vof_c; %save mean VOF
    else %user wants to load a VOF from file instead
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
        [file, path] = uigetfile({'*.*'},'Select VOF for analysis.'); %invoke file select dialog with this title
        if file == 0 %if user cancels dialog
            wrap_text = 'No VOF selected. Analysis terminated.'; %report and terminate
            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
            return
        else %if user actually selects a file
            path_to_vof = fullfile(path, file); %generate full file path
            mean_vof_c = read_this_file(path_to_vof, 1); %call file reader, specify output must be 1-dimensional; function raises error otherwise
            if handles.display_vof %user wants to display the VOF
                vof_axes = handles.axes2;
                if ~handles.display_aif && ~handles.make_plots; vof_axes = handles.axes1; end
                t = 0:handles.tr:(handles.img_size(4)-1)*handles.tr;  %generate time vector for plotting
                plot(t, mean_vof_c, 'k-', 'Parent',vof_axes); xlabel(vof_axes, 't [s]'); %plot
                ylabel(vof_axes, '{C_v}(t)'); title(vof_axes, 'VOF')
            end
        end
    end
    handles.mean_vof_c = mean_vof_c; %save the loaded VOF
    if ~exist('vof_area', 'var'); t = 0:handles.tr:(handles.img_size(4)-1)*handles.tr; %if area under VOF has not been calculated, generate time vector
        vof_area = trapz(t, mean_vof_c); end %and calculate area under VOF
    handles.vof_area = vof_area; %save area
    handles.pve_correction = handles.aif_area/handles.vof_area; %PVE correction factor is AUC(AIF) / AUC(VOF)
end
guidata(hObject, handles) %update app data

%% FIND ARRIVAL TIME OF AIF
if (0)% Redundant
    aif_arrival_time = find_aif_arrival_time(handles);
end

%% DISPLAY DSC_DATA AND REQUEST USER TO SELECT REGION TO ANALYSE
%Invoke question dialog to ask user if user wants to define ROI to restrict
%analysis to particular regions of the brain.
opts.Interpreter = 'tex';
opts.Default = 'No';
response = questdlg('\fontsize{10}Program allows ROI placement to restrict analysis to user-defined regions. Continue with ROI placement?', ...
    'DECONVOLVER', opts);
switch response
    case 'No'; wants_roi = false;
    case 'Yes'; wants_roi = true;
    case 'Cancel'; set(handles.edit7, 'String', 'Analysis terminated by user. Re-run to select data.', 'ForeGroundColor', 'r'); return
    case ''; set(handles.edit7, 'String', 'Analysis terminated by user. Rerun to select data.', 'ForeGroundColor', 'r'); return
end
%__________________________________________________________________________
if wants_roi
    %Turn on image-related objects on axes1; off on axes2;
    set(handles.slider1, 'Visible', 'on')%slice change on axes1
    set(handles.slider2, 'Visible', 'off')%slice change on axes2
    set(handles.slider6, 'Visible', 'on') %upper window setting on axes1
    set(handles.slider10, 'Visible', 'on')%lower window setting on axes1
    set(handles.slider7, 'Visible', 'off')%upper window setting on axes2
    set(handles.slider11, 'Visible', 'off')%lower window setting on axes2
    set(handles.slider8, 'Visible', 'on')%horizontal scrolling on axes1
    set(handles.slider9, 'Visible', 'off')%horizontal scrolling on axes2
    set(handles.pushbutton15, 'Visible', 'on')%colormap on axes1
    set(handles.pushbutton16, 'Visible', 'off')%colormap on axes2
    set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
    set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes2
    set(handles.togglebutton2, 'Visible', 'on')%colorbar on axes1
    set(handles.togglebutton3, 'Visible', 'off')%colorbar on axes2
    %__________________________________________________________________________
    handles.choosing_region_to_analyse = true; %alert objects user is now selecting region to analyse
    initial_slice = ceil(handles.img_size(3)/2); %show middle slice as initial slice
    if ~exist('t_min_signal', 'var') %if the time-point for minimum signal does not exist
        current_state = handles.make_plots;
        handles.make_plots = false;
        [~, t_min_signal] = whole_brain_curve_gui(handles); %extract it by generating whole-brain-curve
        handles.t_min_signal = t_min_signal;
        handles.make_plots = current_state;
    else
        t_min_signal = handles.t_min_signal; %save time point for minimum signal
    end
    if isfield(handles, 'mask') %if mask has been generated
        handles.dsc_data_filtered = handles.dsc_data.*0;
        for sl = 1:handles.img_size(3) %filter image to be displayed using the mask
            for tp = 1:handles.img_size(4)
                handles.dsc_data_filtered(:,:,sl, tp) = squeeze(handles.dsc_data(:,:,sl, tp)).*squeeze(handles.mask(:,:,sl));
            end
        end
    end
    initial_slice_image = squeeze(handles.dsc_data_filtered(:,:,initial_slice,t_min_signal));
    clim_im = [0   max(initial_slice_image(isfinite(initial_slice_image)))]; %color limit, can be changed using sliders
    handles.dsc_data_image = imshow(zeros(handles.img_size(1), handles.img_size(2)),clim_im, 'colormap', gray, 'Parent', handles.axes1);
    title(handles.axes1, ['DSC signal: Slice ' num2str(initial_slice), ' : Time point ' num2str(t_min_signal)])
    set(handles.dsc_data_image, 'CData', initial_slice_image ) %display initial slice
    set(handles.slider1, 'Min', 1, 'Max', handles.img_size(3), ... %setup slider for changing slice
        'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', initial_slice)
    set(handles.slider8, 'Visible', 'on')%horizontal scrolling on axes1
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
    guidata(hObject, handles)
    %__________________________________________________________________________
    n_rois = 128*128; %maximum number of ROIs
    no_roi = false; %will be true if user cancels ROI placement
    this_pos = 0; %keeps track of how many ROIs have been drawn
    roi_type = 'Freehand'; %choose initially a freehand ROI
    for this_roi = 1:n_rois
        if this_roi > 1 %if user has already drawn a ROI
            wrap_text = ['Click to add another ROI. Press: Shift + [F (Freehand), R (rectangle), C (circle), E (ellipse)]. ',...
                'Press Esc to terminate ROI selection. Press Enter to continue with current selection(s).']; %display this message
            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
            enough = waitforbuttonpress; %wait for key press or mouse click
            if enough %if user presses key
                this_key = double(get(handles.figure1, 'CurrentCharacter')); %check what key it was
                switch this_key
                    case 27;  no_roi = true;  %Esc; user wants to skip ROI placement
                    case 82; roi_type = 'Rectangle'; %Shift R
                    case 67; roi_type = 'Circle'; %Shift C
                    case 69; roi_type = 'Ellipse'; %Shift E
                    case 70; roi_type = 'Freehand'; %Shift F
                    case 13; break; %Enter
                end
            end
        end
        if ~no_roi %if user has not chosen to skip ROI placement
            wrap_text = ['Begin drawing of ROI.  Double-click in ROI when done. '...
                'Use slider to change displayed slice.' ]; %display this message
            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
            blink_this(handles.edit7, 'b') %make the text blink
            switch roi_type %call the correct ROI drawing function depending on the current selection
                case 'Freehand'; dsc_roi = drawfreehand(handles.axes1, 'Color', 'yellow', 'LineWidth', 0.1);
                case 'Rectangle'; dsc_roi = drawrectangle(handles.axes1, 'Color', 'yellow', 'LineWidth', 0.1);
                case 'Circle'; dsc_roi = drawcircle(handles.axes1, 'Color', 'yellow', 'LineWidth', 0.1);
                case 'Ellipse'; dsc_roi = drawellipse(handles.axes1, 'Color', 'yellow', 'LineWidth', 0.1);
            end
            dsc_roi_pos = customWait(dsc_roi); %wait for double-click in drawn ROI
            if ishandle(dsc_roi) %if double-clicked ROI still exists
                set(dsc_roi, 'Color', 'green') %colour it green
                try set(dsc_roi, 'FaceSelectable', 0); catch; end %try to make that ROI undeletable
            end
            this_pos = this_pos + 1;%increase ROI counter
            handles.dsc_rois.(['roi_' num2str(this_pos)]) = dsc_roi; %save ROI to this structure
        else %if user has chosen to skip ROI placement
            drawn_rois = findobj(handles.axes1, 'Type', 'images.roi'); %find all drawn ROIS
            delete(drawn_rois) %and delete them
            break %exit loop
        end
    end
    
    if ~no_roi %if user did not skip ROI placement
        if ~isfield(handles, 'mask'); handles.mask = true(handles.img_size(1:3)); end     %if user chose not to generate mask, create a mask true everywhere
        mask_roi = false(size(handles.mask)); %initialise ROI mask, false everywhere
        pb = waitbar(0, 'Incorporating user-defined region(s)...', 'Name', 'DECONVOLVER'); %create waitbar; it will take a while to incorporate selected pixels
        progress_step = 0; %initial progress bar excursion
        tot_length = handles.img_size(1)*handles.img_size(2); %total progress bar excursion
        for x = 1:handles.img_size(1)
            for y = 1:handles.img_size(2) %for each pixel
                pixel_is_included = false; %initially assume no roi membership
                for a_roi = 1:this_pos %go through all selected rois
                    if ishandle(handles.dsc_rois.(['roi_' num2str(a_roi)])) && inROI(handles.dsc_rois.(['roi_' num2str(a_roi)]), x, y) %if pixel exists in at least one of those rois
                        pixel_is_included = true; %set roi membership to true
                    end
                end
                if pixel_is_included %if pixel is member of some roi
                    mask_roi(y,x,:) = handles.mask(y,x,:)*true; %mask value for that pixel is true iff original mask value was also true
                else
                    mask_roi(y,x,:) = false; %otherwise mask value is false
                end
                progress_step = progress_step + 1; %increase progress
                waitbar(progress_step/tot_length, pb) %update progress bar
            end
        end
        delete(pb) %delete progress bar
    end
else %if user does not want ROI
    no_roi = true;
end
if no_roi; handles.mask = mask; else; handles.mask = mask_roi; end %update mask
if ~handles.axes1_vaccant; handles.choosing_region_to_analyse = false; end %if axes1 has not been reserved, keep DSC image and selected ROIs on that axes.
guidata(hObject, handles) %update app data
%% BEZIER CURVE DECONVOLUTION
pause(0.5)
if handles.BzD %if user chose Bezier Curve deconvolution
    if handles.edit_algo_setup %and chose to edit kappa and deconvolution slice range
        handles = get_deconvolution_settings(handles, 'bzd'); %get these settings
        if handles.user_cancelled; handles.user_cancelled = []; return; end %if user cancelled input of these settings, terminate execution
        guidata(hObject, handles) %update app data
    end
    
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
    if handles.live_display %if user wants to display results during analysis
        handles = ready_the_axes(handles); %prepare the axes for this
        params = bezier_deconvolve_gui_vis(handles); %perform BzD with live data visualisation
    else
        params = bezier_deconvolve_gui(handles); %otherwise without live data visualisation
    end
    
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
    if handles.edit_algo_setup %and wants to edit kappa, OI, Psvd, etc
        if handles.do_oSVD; caller = 'osvd'; end %let 'get_deconvolution_settings' which algorithm has been chosen
        if handles.do_cSVD; caller = 'csvd'; end
        if handles.do_sSVD; caller = 'ssvd'; end
        handles = get_deconvolution_settings(handles, caller); %get the relevant settings
        if handles.user_cancelled; handles.user_cancelled = false; return; end %if user cancelled input dialog, return
        guidata(hObject, handles)
    end
    pause(0.5)
    %create progress bar with correct text depending on chosen algorithm____
    if handles.do_oSVD; wrap_text = 'Running oSVD deconvolution ...'; end
    if handles.do_cSVD; wrap_text = 'Running cSVD deconvolution ...'; end
    if handles.do_sSVD; wrap_text = 'Running sSVD deconvolution ...';  end
    %_________________________________________________________________________
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b') %display same message on main window
    pause(0.5)
    %     handles.p_bar = waitbar(0, [p_bar_id '0% complete'], 'Name', 'DECONVOLVER');
    
    childlist = get(handles.figure1, 'Children');
    for c = 1:numel(childlist)
        try set(childlist(c), 'Enable', 'off'); catch; end
    end
    %         set(handles.figure1, 'pointer', 'watch')
    
    pause(0.5)
    if handles.live_display %user wants to diaply results during analysis
        handles = ready_the_axes(handles); %prepare the axes
        pause(0.5)
        guidata(hObject, handles)
        [fitd_r_svd, fitd_cbf_svd, fitd_delay_svd] = SVD_deconvolution_gui_vis(handles); %call the function that displays live results
    else
        if handles.do_oSVD
            [fitd_r_svd, fitd_cbf_svd, fitd_delay_svd] = oSVD_deconvolution_gui(handles); %otherwise call the conventional function
        end
        if handles.do_cSVD
            [fitd_r_svd, fitd_cbf_svd, fitd_delay_svd] = cSVD_deconvolution_gui(handles); %otherwise call the conventional function
        end
        if handles.do_sSVD
            [fitd_r_svd, fitd_cbf_svd, fitd_delay_svd] = sSVD_deconvolution_gui(handles); %otherwise call the conventional function
        end
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
if handles.calculate_oef || handles.calculate_cmro2
    waitbar(6/7, handles.p_bar, 'Calculating OEF and/or CMRO2')
    [handles.oef_BzD, handles.cmro2_BzD, handles.oef_svd, handles.cmro2_svd] = oef_cmro2_gui(handles);
end

%Done______________________________________________________________________
waitbar(1, handles.p_bar, 'Done.')
delete(handles.p_bar)
set(handles.edit7, 'String', 'Analysis is complete.' ,'ForeGroundColor', 'blue')

%% DISPLAY RESULTS
if handles.post_display
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
    
    if handles.show_delay
        if handles.BzD; handles.delay_data =  handles.fitd_delay; else; handles.delay_data =  handles.fitd_delay_svd; end
        set(handles.delay_image, 'CData', imrotate(handles.delay_data(:,:,1), 90))
        title(handles.delay_axes, 'Delay [s] : Slice number 1')
        set(handles.delay_slider, 'Min', 1, 'Max', handles.img_size(3), ...
            'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', 1)
        clim_max = double(max(handles.delay_data(isfinite(handles.delay_data)))) + 1E-10;
        if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/999; end
        if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end
        set(handles.delay_window_slider_up, 'Visible', 'on') %color setting
        set(handles.delay_window_slider_up, 'Min', 10E-10, 'Max', clim_max, ...
            'SliderStep', [step_small, step_big], 'Value', clim_max)
        set(handles.delay_window_slider_low, 'Visible', 'on') %color setting
        set(handles.delay_window_slider_low, 'Min', 0, 'Max', clim_max-1E-10, ...
            'SliderStep', [step_small, step_big], 'Value', 0)
        guidata(hObject, handles)
    end
    
    if handles.show_oef
        if handles.BzD; handles.oef_data =  handles.oef_BzD; else; handles.oef_data =  handles.oef_svd; end
        set(handles.oef_image, 'CData', imrotate(handles.oef_data(:,:,1), 90))
        title(handles.oef_axes, 'OEF : Slice number 1')
        set(handles.oef_slider, 'Min', 1, 'Max', handles.img_size(3), ...
            'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', 1)
        clim_max = double(max(handles.oef_data(isfinite(handles.oef_data)))) + 1E-10;
        if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/999; end
        if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end
        set(handles.oef_window_slider_up, 'Visible', 'on') %color setting
        set(handles.oef_window_slider_up, 'Min', 10E-10, 'Max', clim_max, ...
            'SliderStep', [step_small, step_big], 'Value', clim_max)
        set(handles.oef_window_slider_low, 'Visible', 'on') %color setting
        set(handles.oef_window_slider_low, 'Min', 0, 'Max', clim_max-1E-10, ...
            'SliderStep', [step_small, step_big], 'Value', 0)
        guidata(hObject, handles)
    end
    if handles.show_cmro2
        if handles.BzD; handles.cmro2_data =  handles.cmro2_BzD; else; handles.cmro2_data =  handles.cmro2_svd; end
        set(handles.cmro2_image, 'CData', imrotate(handles.cmro2_data(:,:,1),90))
        title(handles.oef_axes, 'CMRO2 [ml/100g/min] : Slice number 1')
        set(handles.cmro2_slider, 'Min', 1, 'Max', handles.img_size(3), ...
            'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', 1)
        clim_max = double(max(handles.cmro2_data(isfinite(handles.cmro2_data)))) + 1E-10;
        if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/999; end
        if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end
        set(handles.cmro2_window_slider_up, 'Visible', 'on') %color setting
        set(handles.cmro2_window_slider_up, 'Min', 10E-10, 'Max', clim_max, ...
            'SliderStep', [step_small, step_big], 'Value', clim_max)
        set(handles.cmro2_window_slider_low, 'Visible', 'on') %color setting
        set(handles.cmro2_window_slider_low, 'Min', 0, 'Max', clim_max-1E-10, ...
            'SliderStep', [step_small, step_big], 'Value', 0)
        guidata(hObject, handles)
    end
    if handles.show_r10
        if handles.BzD; handles.r10_data =  handles.r10_BzD; else; handles.r10_data =  handles.r10_svd; end
        set(handles.r10_image, 'CData', imrotate(handles.r10_data(:,:,1),90))
        title(handles.r10_axes, 'R10 [s] : Slice number 1')
        set(handles.r10_slider, 'Min', 1, 'Max', handles.img_size(3), ...
            'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', 1)
        clim_max = double(max(handles.r10_data(isfinite(handles.r10_data)))) + 1E-10;
        if clim_max <= 99; step_small = clim_max/100; step_big = clim_max/99; end
        if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end
        set(handles.r10_window_slider_up, 'Visible', 'on') %color setting
        set(handles.r10_window_slider_up, 'Min', 10E-10, 'Max', clim_max, ...
            'SliderStep', [step_small, step_big], 'Value', clim_max)
        set(handles.r10_window_slider_low, 'Visible', 'on') %color setting
        set(handles.r10_window_slider_low, 'Min', 0, 'Max', clim_max-1E-10, ...
            'SliderStep', [step_small, step_big], 'Value', 0)
        guidata(hObject, handles)
    end
    if handles.show_r50
        if handles.BzD; handles.r50_data =  handles.r50_BzD; else; handles.r50_data =  handles.r50_svd; end
        set(handles.r50_image, 'CData', imrotate(handles.r50_data(:,:,1),90))
        title(handles.r50_axes, 'R50 [s] : Slice number 1')
        set(handles.r50_slider, 'Min', 1, 'Max', handles.img_size(3), ...
            'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', 1)
        clim_max = double(max(handles.r50_data(isfinite(handles.r50_data)))) + 1E-10;
        if clim_max <= 99; step_small = clim_max/100; step_big = clim_max/99; end
        if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end %$$$$$$$
        set(handles.r50_window_slider_up, 'Visible', 'on') %color setting
        set(handles.r50_window_slider_up, 'Min', 10E-10, 'Max', clim_max, ...
            'SliderStep', [step_small, step_big], 'Value', clim_max)
        set(handles.r50_window_slider_low, 'Visible', 'on') %color setting
        set(handles.r50_window_slider_low, 'Min', 0, 'Max', clim_max-1E-10, ...
            'SliderStep', [step_small, step_big], 'Value', 0)
        guidata(hObject, handles)
    end
    if handles.plot_residue_funcs
        handles.r_slice = handles.slice_range(5);
        title(handles.rt_axes, ['Residue functions : Slice number ' num2str(handles.r_slice)])
        guidata(hObject, handles)
        plot_residue_functions_gui(handles)
        set(handles.rt_slider, 'Min', 1, 'Max', handles.img_size(3), ...
            'SliderStep', [1, 1]/(handles.img_size(3) - 1), 'Value', 1)
        guidata(hObject, handles)
    end
end


%__________________________notify__________________________________________
if handles.notify_when_done
    for bp = 1:3; beep; pause(0.5); beep; end
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
function slider1_Callback(hObject, eventdata, handles)
index = round(get(hObject, 'Value'));
if handles.data_onto_axes1 %user has loaded data onto axes1
    if length(handles.size_axes1_data) == 4
        if handles.plotting_4d_data
            title(handles.axes1, [handles.axes1_filename ': Slice ' num2str(index)], 'Interpreter', 'none')
            plot_4d_data(handles, handles.axes1_data, index, handles.axes1, handles.axes1_filename)
        else
            time_point = round(get(handles.slider8,'Value'));
            this_slice = round(get(handles.slider1,'Value'));
            set(handles.axes1_image, 'CData', squeeze(handles.axes1_data(:,:,this_slice, time_point))  )
            title(handles.axes1, [handles.axes1_filename ': Slice ' num2str(this_slice) ' : Time-point ' num2str(time_point)],...
                'Interpreter', 'None')
        end
    else
        set(handles.axes1_image, 'CData', squeeze(handles.axes1_data(:,:,index))  )
        title(handles.axes1, [handles.axes1_filename ': Slice ' num2str(index)], 'Interpreter', 'None')
    end
else %no data loaded onto axes1
    
    if handles.choosing_region_to_analyse
        time_indx = round(get(handles.slider8, 'Value'));
        set(handles.dsc_data_image, 'CData', squeeze(handles.dsc_data_filtered(:,:,index,time_indx))  )
        title(handles.axes1, ['DSC signal: Slice ' num2str(index), ' : Time point ', num2str(time_indx)])
    else
        
        if handles.data_is_for_aif
            aif_time_point = round(get(handles.slider8, 'Value'));
            set(handles.aif_conc_image, 'CData', squeeze(handles.aif_conc_data(:,:,index,aif_time_point)))
            title(handles.axes1, ['Concentration: Slice ' num2str(index) ' : Time point ' num2str(aif_time_point)])
        else
            if handles.data_is_for_vof
                vof_time_point = round(get(handles.slider8, 'Value'));
                set(handles.vof_conc_image, 'CData', squeeze(handles.vof_conc_data(:,:,index,vof_time_point)))
                title(handles.axes1, ['Concentration: Slice ' num2str(index) ' : Time point ' num2str(vof_time_point)])
            else
                
                
                if handles.show_cbf && handles.cbf_slider == handles.slider1
                    set(handles.cbf_image, 'CData', imrotate(handles.cbf_data(:,:,index), 90) );
                    title(handles.axes1, ['CBF [ml/100g/min] : Slice number ' num2str(index)])
                end
                if handles.show_mtt && handles.mtt_slider == handles.slider1
                    set(handles.mtt_image, 'CData', imrotate(handles.mtt_data(:,:,index), 90) );
                    title(handles.axes1, ['MTT [s] : Slice number ' num2str(index)])
                end
                if handles.show_delay && handles.delay_slider == handles.slider1
                    set(handles.delay_image, 'CData', imrotate(handles.delay_data(:,:,index), 90));
                    title(handles.axes1, ['Delay [s] : Slice number ' num2str(index)])
                end
                if handles.plot_residue_funcs && handles.rt_slider == handles.slider1
                    handles.r_slice = index;
                    title(handles.axes1, ['Residue functions: Slice number ' num2str(index)])
                    plot_residue_functions_gui(handles)
                end
                if handles.show_oef && handles.oef_slider == handles.slider1
                    set(handles.oef_image, 'CData', imrotate(handles.oef_data(:,:,index), 90));
                    title(handles.axes1, ['OEF : Slice number ' num2str(index)])
                end
                if handles.show_cmro2 && handles.cmro2_slider == handles.slider1
                    set(handles.cmro2_image, 'CData', imrotate(handles.cmro2_data(:,:,index), 90));
                    title(handles.axes1, ['CMRO2 [ml/100g/min] : Slice number ' num2str(index)])
                end
                if handles.show_r10 && handles.r10_slider == handles.slider1
                    set(handles.r10_image, 'CData', imrotate(handles.r10_data(:,:,index), 90));
                    title(handles.axes1, ['R10 [s] : Slice number ' num2str(index)])
                end
                if handles.show_r50 && handles.r50_slider == handles.slider1
                    set(handles.r50_image, 'CData', imrotate(handles.r50_data(:,:,index), 90));
                    title(handles.axes1, ['R50 [s] : Slice number ' num2str(index)])
                end
            end
        end
    end
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

if handles.data_onto_axes2 %user has loaded data onto axes1
    if length(handles.size_axes2_data) == 4
        if handles.plotting_4d_data
            title(handles.axes2, [handles.axes2_filename ': Slice ' num2str(index)], 'Interpreter', 'none')
            plot_4d_data(handles, handles.axes2_data, index, handles.axes2, handles.axes2_filename)
        else
            time_point = round(get(handles.slider9,'Value'));
            this_slice = round(get(handles.slider2,'Value'));
            set(handles.axes2_image, 'CData', squeeze(handles.axes2_data(:,:,this_slice, time_point))  )
            title(handles.axes2, [handles.axes2_filename ': Slice ' num2str(this_slice) ' : Time-point ' num2str(time_point)],...
                'Interpreter', 'None')
        end
    else
        set(handles.axes2_image, 'CData', squeeze(handles.axes2_data(:,:,index))  )
        title(handles.axes2, [handles.axes2_filename ': Slice ' num2str(index)], 'Interpreter', 'None')
    end
else
    
    if handles.data_is_for_aif
        aif_time_point = handles.aif_handles.aif_time_point;
        set(handles.aif_conc_image, 'CData', squeeze(handles.aif_conc_data(:,:,index,aif_time_point)))
        title(handles.axes2, ['Concentration: Slice number ' num2str(index)])
    else
        if handles.data_is_for_vof
            vof_time_point = handles.vof_handles.vof_time_point;
            set(handles.vof_conc_image, 'CData', squeeze(handles.vof_conc_data(:,:,index,vof_time_point)))
            title(handles.axes2, ['Concentration: Slice number ' num2str(index)])
        else
            
            if handles.show_cbf && handles.cbf_slider == handles.slider2
                set(handles.cbf_image, 'CData', imrotate(handles.cbf_data(:,:,index), 90) );
                title(handles.axes2, ['CBF [ml/100g/min] : Slice number ' num2str(index)])
            end
            if handles.show_mtt && handles.mtt_slider == handles.slider2
                set(handles.mtt_image, 'CData', imrotate(handles.mtt_data(:,:,index), 90) );
                title(handles.axes2, ['MTT [s] : Slice number ' num2str(index)])
            end
            if handles.show_delay && handles.delay_slider == handles.slider2
                set(handles.delay_image, 'CData', imrotate(handles.delay_data(:,:,index), 90));
                title(handles.axes2, ['Delay [s] : Slice number ' num2str(index)])
            end
            if handles.plot_residue_funcs && handles.rt_slider == handles.slider2
                handles.r_slice = index;
                title(handles.axes2, ['Residue functions: Slice number ' num2str(index)])
                plot_residue_functions_gui(handles)
            end
            if handles.show_oef && handles.oef_slider == handles.slider2
                set(handles.oef_image, 'CData', imrotate(handles.oef_data(:,:,index), 90));
                title(handles.axes2, ['OEF : Slice number ' num2str(index)])
            end
            if handles.show_cmro2 && handles.cmro2_slider == handles.slider2
                set(handles.cmro2_image, 'CData', imrotate(handles.cmro2_data(:,:,index), 90));
                title(handles.axes2, ['CMRO2 [ml/100g/min] : Slice number ' num2str(index)])
            end
            if handles.show_r10 && handles.r10_slider == handles.slider2
                set(handles.r10_image, 'CData', imrotate(handles.r10_data(:,:,index), 90));
                title(handles.axes2, ['R10 [s] : Slice number ' num2str(index)])
            end
            if handles.show_r50 && handles.r50_slider == handles.slider2
                set(handles.r50_image, 'CData', imrotate(handles.r50_data(:,:,index), 90));
                title(handles.axes2, ['R50 [s] : Slice number ' num2str(index)])
            end
        end
    end
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in DSC-DATA popuo menu
function popupmenu6_Callback(hObject, eventdata, handles)
data_choices = cellstr(get(hObject,'String'));
data_choice = data_choices{get(hObject,'Value')};
handles.load_existing_data = false;
handles.simulate_data = false;
switch data_choice
    case 'Load'
        handles.load_existing_data = true;
    case 'Simulate'
        handles.simulate_data = true;
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in EDIT SETUP FOR DSC-MRI DATA.
function checkbox24_Callback(hObject, eventdata, handles)
handles.edit_data_setup = get(hObject,'Value');
guidata(hObject, handles)


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


% --- Executes on button press in EDIT ALGORITHM SETUP.
function checkbox25_Callback(hObject, eventdata, handles)
handles.edit_algo_setup = get(hObject,'Value');
guidata(hObject, handles)


% --- Executes on selection change in DISPLAY RESULTS: WHEN.
function popupmenu8_Callback(hObject, eventdata, handles)
display_choices = cellstr(get(hObject,'String'));
display_choice = display_choices{get(hObject,'Value')};
switch display_choice
    case 'After analysis'
        handles.post_display = true;
        handles.live_display = false;
    case 'During analysis'
        handles.live_display = true;
        handles.post_display = true;
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in SAVE RESULTS: FORMAT.
function popupmenu10_Callback(hObject, eventdata, handles)
format_choices = cellstr(get(hObject,'String'));
format_choice = format_choices{get(hObject,'Value')};
handles.save_format = format_choice;
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
% --- Executes on button press in LOAD DATA ONTO AXES1.
function pushbutton8_Callback(hObject, eventdata, handles)
cla(handles.axes1, 'reset')
set(handles.axes1,'XTick', [], 'YTick', [])
set(handles.axes1, 'box', 'on')
handles.plotting_4d_data = false;
set(handles.slider1, 'Visible', 'on')
[file, path] = uigetfile('*');
if file == 0
else
    handles.data_onto_axes1 = true;
    path_to_axes1_data = fullfile(path, file);
    wrap_text = ['Visualising data from ' path_to_axes1_data];
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    [~, file_name, file_ext] = fileparts(path_to_axes1_data);
    %recognise files generated by DECONVOLVER and read them accordingly
    switch file_name
        case {'cbf_BzD', 'cbf_SVD', 'cbv_as_cbf_times_mtt_BzD', 'cbv_as_cbf_times_mtt_SVD', 'cbv_without_deconvolution', 'cmro2_BzD', 'cmro2_SVD',...
                'delay_BzD', 'delay_SVD', 'mtt_as_area_under_r_BzD', 'mtt_as_area_under_r_SVD', 'mtt_as_cbv_over_cbf_BzD', 'mtt_as_cbv_over_cbf_SVD',...
                'oef_BzD', 'oef_SVD', 'ttp_from_calculated_signal_SVD', 'ttp_from_fitted_signal_BzD', 'ttp_without_deconvolution'}; num_dims = 3;
        case {'residue_functions_BzD', 'residue_functions_SVD', 'control_points_BzD', 'disp_kernel_BzD'}; num_dims = 4;
        otherwise num_dims = false;
    end
    axes1_pre_data = read_this_file(path_to_axes1_data, num_dims);
    if isstruct(axes1_pre_data)
        axes1_data = squeeze(axes1_pre_data.the_data);
        handles.header_info = axes1_pre_data.the_info;
    else
        if isequal(axes1_pre_data, false)
            return
        else
            axes1_data = squeeze(axes1_pre_data); %remove singular dimensions as they will confuse sliders
        end
    end
    size_axes1_data = size(axes1_data);
    handles.size_axes1_data = size_axes1_data;
    handles.axes1_data = axes1_data;
    handles.path_to_axes1_data = path_to_axes1_data;
    handles.axes1_filename = file_name;
    
    clim_max = double(max(axes1_data(:)))+1E-10;
    if clim_max <= 99; step_small = clim_max/100; step_big = clim_max/99; end
    if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end
    if length(size_axes1_data) == 2 %2D or 1D
        if ~isempty(size_axes1_data(size_axes1_data == 1)) %1D data
            set(handles.slider1, 'Visible', 'off')
            set(handles.slider6, 'Visible', 'off')
            set(handles.slider10, 'Visible', 'off')
            set(handles.slider8, 'Visible', 'off')
            set(handles.pushbutton15, 'Visible', 'off')
            set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
            set(handles.togglebutton2, 'Visible', 'off')
            opts.Interpreter = 'tex';
            prompt = {['\color{blue} \fontsize{10}You have loaded 1D data (',...
                num2str(size_axes1_data(1)), '\times', num2str(size_axes1_data(2)),') ',...
                'onto axes1. Program will display  your data as a plot in the Cartesian plane.', ...
                sprintf(['\n' '\n']), '\color{blue} \fontsize{10}If applicable, provide x-axis sampling interval below.'...
                ' Leave edit field unchanged otherwise.']};
            defaults = {'1'};
            dims = [1 100];
            dlgtitle = 'Visualisation of 1D data';
            response = inputdlg(prompt, dlgtitle, dims, defaults, opts);
            if ~isempty(response) %Clicked OK
                d_x = str2double(response{1});
                x = 0:d_x:d_x*(max(size_axes1_data)-1);
                plot(x, axes1_data, 'k-','LineWidth',1, 'Parent', handles.axes1)
                title(handles.axes1, file_name, 'Interpreter', 'none')
            end
        else %2D data
            set(handles.slider1, 'Visible', 'off')
            set(handles.slider6, 'Visible', 'on')
            set(handles.slider10, 'Visible', 'on')
            set(handles.slider8, 'Visible', 'off')
            set(handles.pushbutton15, 'Visible', 'on')
            set(handles.pushbutton17, 'Visible', 'on')%rotate image on axes1
            set(handles.togglebutton2, 'Visible', 'on')
            handles.axes1_image = imshow(axes1_data, 'colormap', gray, 'Parent', handles.axes1);
            handles.axes1.CLim = [0 clim_max];
            title(handles.axes1, file_name, 'Interpreter', 'none')
            set(handles.slider1, 'Visible', 'off')
            set(handles.slider6, 'Visible', 'on') %color setting
            set(handles.slider6, 'Min', 10E-10, 'Max', clim_max, ...
                'SliderStep', [step_small, step_big], 'Value', clim_max)
            set(handles.slider10, 'Visible', 'on') %color setting
            set(handles.slider10, 'Min', 0, 'Max', clim_max-1E-10, ...
                'SliderStep', [step_small, step_big], 'Value', 0)
        end
    end
    
    if length(size_axes1_data) == 3 %3D data
        if ( strcmp(file_ext, '.png') || strcmp(file_ext, '.PNG') ...
                || strcmp(file_ext, '.jpg') || strcmp(file_ext, '.JPG') || ...
                strcmp(file_ext, '.JPEG') || strcmp(file_ext, '.gif') || strcmp(file_ext, '.ico'))
            handles.axes1_image = imshow(axes1_data, 'Parent', handles.axes1); %ensure correct visualisation of rgb photos
            set(handles.slider1, 'Visible', 'off')
            set(handles.slider6, 'Visible', 'on')
            set(handles.slider10, 'Visible', 'on')
            set(handles.slider8, 'Visible', 'off')
            set(handles.pushbutton15, 'Visible', 'on')
            set(handles.pushbutton17, 'Visible', 'on')%rotate image on axes1
            set(handles.togglebutton2, 'Visible', 'on')
            title(handles.axes1, file_name, 'Interpreter', 'none')
            guidata(hObject, handles)
        else
            set(handles.slider1, 'Visible', 'on')
            set(handles.slider6, 'Visible', 'on')
            set(handles.slider10, 'Visible', 'on')
            set(handles.slider8, 'Visible', 'off')
            set(handles.pushbutton15, 'Visible', 'on')
            set(handles.pushbutton17, 'Visible', 'on')%rotate image on axes1
            set(handles.togglebutton2, 'Visible', 'on')
            handles.axes1_image = imshow(squeeze(axes1_data(:,:,1)), 'colormap', gray, 'Parent', handles.axes1);
            handles.axes1.CLim = [0 clim_max];
            set(handles.slider1, 'Min', 1, 'Max', size_axes1_data(3), ...
                'SliderStep', [1, 1]/(size_axes1_data(3) - 1), 'Value', 1)
            set(handles.slider6, 'Visible', 'on') %color setting
            set(handles.slider6, 'Min', 10E-10, 'Max', clim_max, ...
                'SliderStep', [step_small, step_big], 'Value', clim_max)
            set(handles.slider10, 'Visible', 'on') %color setting
            set(handles.slider10, 'Min', 0, 'Max', clim_max-1E-10, ...
                'SliderStep', [step_small, step_big], 'Value', 0)
            title(handles.axes1, [file_name, ': Slice 1'], 'Interpreter', 'none')
        end
    end
    
    if length(size_axes1_data) == 4 %nasty
        opts.Interpreter = 'tex';
        prompt = {['\color{blue} \fontsize{10}You have loaded 4D data (',...
            num2str(size_axes1_data(1)), '\times', num2str(size_axes1_data(2)), '\times', num2str(size_axes1_data(3)), ...
            '\times', num2str(size_axes1_data(4)), ') ',...
            'onto axes1. Program can either display your data as 1D plots in the Cartesian plane or image sequences.', ...
            ' Specifiy desired display (plots or images).'], [sprintf(['\n' '\n']), ...
            '\color{blue} \fontsize{10}If applicable, provide x-axis/temporal sampling interval below.'...
            ' Leave edit field unchanged otherwise.']};
        defaults = {'images', '1'};
        dims = [1 100];
        dlgtitle = 'Visualisation of 4D data';
        response = inputdlg(prompt, dlgtitle, dims, defaults, opts);
        if ~isempty(response) %Clicked OK
            d_t = str2double(response{2});
            t = 0:d_t:d_t*(max(size_axes1_data(4))-1);
            handles.fourth_dim_vector = t;
            if ~strcmp(response{1}, 'plots') %if user selects image display mode
                %                     title(handles.axes1, file_name, 'Interpreter', 'none')
                set(handles.slider1, 'Visible', 'on')
                set(handles.slider6, 'Visible', 'on')
                set(handles.slider10, 'Visible', 'on')
                set(handles.slider8, 'Visible', 'on')
                set(handles.pushbutton15, 'Visible', 'on')
                set(handles.pushbutton17, 'Visible', 'on')%rotate image on axes1
                set(handles.togglebutton2, 'Visible', 'on')
                handles.axes1_image = imshow(squeeze(axes1_data(:,:,1,1)), 'colormap', gray, 'Parent', handles.axes1);
                handles.axes1.CLim = [0 clim_max];
                set(handles.slider1, 'Min', 1, 'Max', size_axes1_data(3), ... %1st adn 2nd dimensions
                    'SliderStep', [1, 1]/(size_axes1_data(3) - 1), 'Value', 1)
                dummy_data = axes1_data; dummy_data(~isfinite(dummy_data)) = []; dummy_data(isnan(dummy_data)) = [];
                clim_max = double(max(dummy_data(:))) +1E-10;
                clear('dummy_data')
                if clim_max <= 99; step_small = clim_max/100; step_big = clim_max/99; end
                if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end
                set(handles.slider6, 'Visible', 'on') %color setting
                set(handles.slider6, 'Min', 10E-10, 'Max', clim_max, ...
                    'SliderStep', [step_small, step_big], 'Value', clim_max)
                set(handles.slider10, 'Visible', 'on') %color setting
                set(handles.slider10, 'Min', 0, 'Max', clim_max-1E-10, ...
                    'SliderStep', [step_small, step_big], 'Value', 0)
                
                set(handles.slider8, 'Min', 1, 'Max', size_axes1_data(4), ... %4th dimension
                    'SliderStep', [1, 1]/(size_axes1_data(4) - 1), 'Value', 1)
                title(handles.axes1, [file_name, ': Slice 1 : Time-point 1'], 'Interpreter', 'none')
            end
            
            if strcmp(response{1}, 'plots') %if user selects plots display mode
                set(handles.slider1, 'Visible', 'off')
                set(handles.slider6, 'Visible', 'off')
                set(handles.slider10, 'Visible', 'off')
                set(handles.slider8, 'Visible', 'off')
                set(handles.pushbutton15, 'Visible', 'off')
                set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
                set(handles.togglebutton2, 'Visible', 'off')
                %                         wrap_text = 'Waiting for user input...';
                %                         set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                opts.Interpreter = 'tex';
                num_plots = size_axes1_data(1)*size_axes1_data(2);
                prompt = {[sprintf('\n'), ...
                    '\color{blue} \fontsize{10}You are about to plot ',...
                    num2str(size_axes1_data(1)), '\times', num2str(size_axes1_data(2)), ' ( =' , num2str(num_plots), ...
                    ') plots on axes1. ',...
                    '\color{red} \fontsize{10}This is highly unrecommended.', ...
                    '\color{blue} \fontsize{10} Program allows selection of one or more ROIs from which to plot',...
                    '\color{blue} \fontsize{10} Press ''Define ROI'' to load image on which to place ROI(s).',...
                    sprintf(['\n' '\n']), ...
                    '\color{magenta} \fontsize{10}WARNING: Image will be loaded onto axes2. Axes2 will be cleared. ',...
                    sprintf('\n'), ...
                    '\color{blue} \fontsize{10}Select action.']};
                opts.Default = 'Define ROI';
                %                         dims = [1 100];
                dlgtitle = 'Visualisation of 4D data';
                response2 = questdlg(prompt, dlgtitle, 'Define ROI','Do not define ROI', 'Cancel', opts);
                switch response2 %user did not cancel dialog
                    case 'Define ROI' %user chose to define ROI
                        cla(handles.axes2, 'reset')
                        set(handles.axes2,'XTick', [], 'YTick', [])
                        set(handles.axes2, 'box', 'on')
                        roi_type = 'Freehand'; %type of ROI
                        [file, path] = uigetfile('*');
                        [~, file_name2, ~] = fileparts(fullfile(path,file));
                        %recognise files generated by DECONVOLVER and read them accordingly
                        switch file_name2
                            case {'cbf_BzD', 'cbf_SVD', 'cbv_as_cbf_times_mtt_BzD', 'cbv_as_cbf_times_mtt_SVD', 'cbv_without_deconvolution', 'cmro2_BzD', 'cmro2_SVD',...
                                    'control_points_BzD', 'delay_BzD', 'delay_SVD', 'mtt_as_area_under_r_BzD', 'mtt_as_area_under_r_SVD', 'mtt_as_cbv_over_cbf_BzD', 'mtt_as_cbv_over_cbf_SVD',...
                                    'oef_BzD', 'oef_SVD', 'ttp_from_calculated_signal_SVD', 'ttp_from_fitted_signal_BzD', 'ttp_without_deconvolution'}; num_dims = 3;
                            case {'residue_functions_BzD', 'residue_functions_SVD'}; num_dims = 4;
                            otherwise; num_dims = false;
                        end
                        if file == 0 %do nothing if no file is selected
                        else
                            axes2_pre_data = read_this_file(fullfile(path, file),num_dims);
                            if isstruct(axes2_pre_data)
                                axes2_data = squeeze(axes2_pre_data.the_data);  %remove singular dimensions as they will confuse sliders
                                handles.header_info = axes2_pre_data.the_info;
                            else
                                if isequal(axes2_pre_data, false)
                                    return
                                else
                                    axes2_data = squeeze(axes2_pre_data);
                                end
                            end
                            size_axes2_data = size(axes2_data);
                            handles.size_axes2_data = size_axes2_data;
                            handles.axes2_data = axes2_data;
                            clim_max_2 = double(max(axes2_data(:))) + 1E-10;
                            if clim_max_2 <= 99; step_small = clim_max_2/100; step_big = clim_max_2/99; end
                            if clim_max_2 > 99; step_small = 1/clim_max_2; step_big = step_small; end
                            [~, file_name_2,file_ext_2] = fileparts(fullfile(path, file));
                            handles.axes2_filename = file_name_2;
                            handles.data_onto_axes2 = true;
                            guidata(hObject, handles)
                            if length(size_axes2_data) == 2 %loaded 2D image
                                set(handles.slider2, 'Visible', 'off')
                                set(handles.slider7, 'Visible', 'on')
                                set(handles.slider11, 'Visible', 'on')
                                set(handles.slider9, 'Visible', 'off')
                                set(handles.pushbutton16, 'Visible', 'on')
                                set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes1
                                set(handles.togglebutton3, 'Visible', 'on')
                                handles.axes1_roi_image = imshow(axes2_data, 'colormap', gray, 'Parent', handles.axes2);
                                handles.axes1.CLim = [0 clim_max_2];
                                set(handles.slider2, 'Visible', 'off') %slice change
                                set(handles.slider7, 'Visible', 'on') %color setting
                                set(handles.slider7, 'Min', 10E-10, 'Max', clim_max_2, ...
                                    'SliderStep', [step_small, step_big], 'Value', clim_max_2)
                                set(handles.slider11, 'Visible', 'on') %color setting
                                set(handles.slider11, 'Min', 0, 'Max', clim_max_2-1E-10, ...
                                    'SliderStep', [step_small, step_big], 'Value', 0)
                                title(handles.axes2, file_name_2, 'Interpreter', 'none')
                                guidata(hObject, handles)
                            end
                            
                            if length(size_axes2_data) == 3 %loaded 3D image
                                if ( strcmp(file_ext_2, '.png') || strcmp(file_ext_2, '.PNG') ...
                                        || strcmp(file_ext_2, '.jpg') || strcmp(file_ext_2, '.JPG') || ...
                                        strcmp(file_ext_2, '.JPEG') || strcmp(file_ext_2, '.gif') || strcmp(file_ext_2, '.ico'))
                                    handles.axes2_image = imshow(axes2_data, 'Parent', handles.axes2); %ensure correct visualisation of rgb photos
                                    set(handles.slider2, 'Visible', 'off')
                                    set(handles.slider7, 'Visible', 'on')
                                    set(handles.slider11, 'Visible', 'on')
                                    set(handles.slider9, 'Visible', 'off')
                                    set(handles.pushbutton16, 'Visible', 'on')
                                    set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes1
                                    set(handles.togglebutton3, 'Visible', 'on')
                                    title(handles.axes2, file_name_2, 'Interpreter', 'none')
                                    guidata(hObject, handles)
                                else
                                    set(handles.slider2, 'Visible', 'on')
                                    set(handles.slider7, 'Visible', 'on')
                                    set(handles.slider11, 'Visible', 'on')
                                    set(handles.slider9, 'Visible', 'off')
                                    set(handles.pushbutton16, 'Visible', 'on')
                                    set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes1
                                    set(handles.togglebutton3, 'Visible', 'on')
                                    handles.axes2_image = imshow(squeeze(axes2_data(:,:,1)), 'colormap', gray, 'Parent', handles.axes2);
                                    handles.axes2.CLim = [0 clim_max_2];
                                    set(handles.slider2, 'Min', 1, 'Max', size_axes2_data(3), ...
                                        'SliderStep', [1, 1]/(size_axes2_data(3) - 1), 'Value', 1)
                                    set(handles.slider7, 'Visible', 'on') %color setting
                                    set(handles.slider7, 'Min', 10E-10, 'Max', clim_max_2, ...
                                        'SliderStep', [step_small, step_big], 'Value', clim_max_2)
                                    set(handles.slider11, 'Visible', 'on') %color setting
                                    set(handles.slider11, 'Min', 0, 'Max', clim_max_2-1E-10, ...
                                        'SliderStep', [step_small, step_big], 'Value', 0)
                                    title(handles.axes2, [file_name_2, ': Slice 1'], 'Interpreter', 'none')
                                    guidata(hObject, handles)
                                end
                            end
                            
                            if length(size_axes2_data) == 4 %loaded 4D image
                                set(handles.slider2, 'Visible', 'on')
                                set(handles.slider7, 'Visible', 'on')
                                set(handles.slider11, 'Visible', 'on')
                                set(handles.slider9, 'Visible', 'on')
                                set(handles.pushbutton16, 'Visible', 'on')
                                set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes1
                                set(handles.togglebutton3, 'Visible', 'on')
                                handles.axes2_image = imshow(squeeze(axes2_data(:,:,1,1)), 'colormap', gray, 'Parent', handles.axes2);
                                handles.axes2.CLim = [0 clim_max_2];
                                set(handles.slider2, 'Min', 1, 'Max', size_axes2_data(3), ... %1st adn 2nd dimensions
                                    'SliderStep', [1, 1]/(size_axes2_data(3) - 1), 'Value', 1)
                                dummy_data = axes2_data; dummy_data(~isfinite(dummy_data)) = []; dummy_data(isnan(dummy_data)) = [];
                                clim_max_2 = double(max(dummy_data(:))) + 1E-10;
                                if clim_max_2 <= 99; step_small = clim_max2/100; step_big = clim_max_2/99; end
                                if clim_max_2 > 99; step_small = 1/clim_max_2; step_big = step_small; end
                                clear('dummy_data')
                                set(handles.slider7, 'Visible', 'on') %color setting
                                set(handles.slider7, 'Min', 10E-10, 'Max', clim_max_2, ...
                                    'SliderStep', [step_small, step_big], 'Value', clim_max_2)
                                set(handles.slider11, 'Visible', 'on') %color setting
                                set(handles.slider11, 'Min', 0, 'Max', clim_max_2-1E-10, ...
                                    'SliderStep', [step_small, step_big], 'Value', 0)
                                set(handles.slider9, 'Min', 1, 'Max', size_axes2_data(4), ... %4th dimension
                                    'SliderStep', [1, 1]/(size_axes2_data(4) - 1), 'Value', 1)
                                title(handles.axes2, [file_name_2, ': Slice 1 : Time-point 1'], 'Interpreter', 'none')
                                guidata(hObject, handles)
                            end
                            %OK, whatever image has been loaded is up
                            %and ready for ROI placement.
                            %----------------------------------------------------------------------------------------------------------------
                            wrap_text = ['Click on the image to begin drawing of freehand ROI(s).  Double-click in ROI when done. '...
                                'Press any key to cancel.'...
                                'Use slider to change displayed slice.' ];
                            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
                            blink_this(handles.edit7, 'b')
                            analyse_all = waitforbuttonpress;
                            n_rois = 128*128;
                            this_pos = 0;
                            handles.plot_colors = {'black', 'green', 'blue', 'red', 'cyan', 'magenta', 'yellow', rgb('Silver'), rgb('DarkGreen'), rgb('SkyBlue'), rgb('DarkRed'), rgb('CadetBlue'),...
                                rgb('Indigo'), rgb('Brown'), rgb('Gray'), rgb('HotPink'), rgb('DarkOrange'), rgb('DodgerBlue'), rgb('RosyBrown'), rgb('Lime'), rgb('OrangeRed'), rgb('GreenYellow'),...
                                rgb('Purple'), rgb('Chocolate'), rgb('SaddleBrown'), rgb('DarkKhaki'), rgb('Pink'), rgb('Salmon'), rgb('Olive'), rgb('LightSteelBlue')};
                            if ~analyse_all %mouse clicked
                                for this_roi = 1:n_rois
                                    if this_roi > 1
                                        wrap_text = ['Click to add another ROI. Press: Shift + [R (rectangle), C (circle), E (ellipse)].',...
                                            'Press Esc to end.', ' Press Spacebar to reset axes.'];
                                        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
                                        enough = waitforbuttonpress;
                                        this_key = double(get(handles.figure1, 'CurrentCharacter'));
                                        if enough %if user presses key
                                            switch this_key
                                                case 27; break %Esc
                                                case 82; roi_type = 'Rectangle'; %Shift R
                                                    %                    case 80; roi_type = 'Point'; %Shift P
                                                case 67; roi_type = 'Circle'; %Shift C
                                                case 69; roi_type = 'Ellipse'; %Shift E
                                                case 70; roi_type = 'Freehand'; %Shift F
                                                case 32 %Spacebar
                                                    current_title = handles.axes2.Title.String;
                                                    cla(handles.axes1, 'reset')
                                                    set(handles.axes1,'XTick', [], 'YTick', [])
                                                    set(handles.axes1, 'box', 'on')
                                                    handles.axes1.Title.String = current_title;
                                                    for rr = 1:this_pos
                                                        if ishandle(handles.fourd_rois.(['roi_' num2str(rr)]))
                                                            delete(handles.fourd_rois.(['roi_' num2str(rr)]))
                                                            curves = fieldnames(handles.(['roi_curves' num2str(rr)]));
                                                            for cc = 1:numel(curves)
                                                                delete(handles.(['roi_curves' num2str(rr)]).(['curve_' num2str(cc)]))
                                                            end
                                                        end
                                                    end
                                                    this_pos = 0;
                                            end
                                        end
                                    end
                                    
                                    wrap_text = ['Select region to analyse. Click on the image to begin drawing of ROI.  Double-click in ROI when done. '...
                                        'Use slider to change displayed slice.' ];
                                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
                                    switch roi_type
                                        case 'Freehand'; fourd_roi = drawfreehand(handles.axes2, 'Color', 'yellow', 'LineWidth', 0.1);
                                        case 'Rectangle'; fourd_roi = drawrectangle(handles.axes2, 'Color', 'yellow', 'LineWidth', 0.1);
                                        case 'Circle'; fourd_roi = drawcircle(handles.axes2, 'Color', 'yellow', 'LineWidth', 0.1);
                                        case 'Ellipse'; fourd_roi = drawellipse(handles.axes2, 'Color', 'yellow', 'LineWidth', 0.1);
                                    end
                                    fourd_roi_pos = customWait(fourd_roi);
                                    this_pos = this_pos + 1;
                                    if ishandle(fourd_roi); set(fourd_roi, 'Color', handles.plot_colors{this_pos}); end
                                    handles.fourd_rois.(['roi_' num2str(this_pos)]) = fourd_roi;
                                    drawn_slice = round(get(handles.slider2, 'Value'));
                                    %if user chose not to generate mask, create a mask true everywhere
                                    curve = 0;
                                    for x = 1:handles.size_axes1_data(1)
                                        for y = 1:handles.size_axes1_data(2) %for each pixel
                                            if inROI(handles.fourd_rois.(['roi_' num2str(this_pos)]), x, y) %if pixel exists in at least one of those rois
                                                curve = curve + 1;
                                                handles.(['roi_curves' num2str(this_pos)]).(['curve_' num2str(curve)]) = plot(t, squeeze(axes1_data(y,x,drawn_slice,:)),...
                                                    'Color', handles.plot_colors{this_pos}, 'Parent', handles.axes1);
                                                hold(handles.axes1, 'on')
                                            end
                                        end
                                    end
                                    title(handles.axes1, file_name, 'Interpreter', 'none')
                                end
                            else %if key pressed, do nothing
                            end
                            %----------------------------------------------------------------------------------------------------------------
                        end
                    case 'Do not define ROI' %user edited roi selection field => does not want to choose ROI
                        handles.this_object = hObject;
                        set(handles.slider1, 'Visible', 'on')
                        set(handles.slider6, 'Visible', 'off')
                        set(handles.slider10, 'Visible', 'off')
                        set(handles.slider8, 'Visible', 'off')
                        set(handles.pushbutton15, 'Visible', 'off')
                        set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
                        set(handles.togglebutton2, 'Visible', 'off')
                        title(handles.axes1, [file_name ': Slice 1'], 'Interpreter', 'none')
                        pause(0.5)
                        plot_4d_data(handles, axes1_data, 1, handles.axes1, file_name)
                        handles.plotting_4d_data = true;
                        set(handles.slider1, 'Min', 1, 'Max', size_axes1_data(3), ... %1st adn 2nd dimensions
                            'SliderStep', [1, 1]/(size_axes1_data(3) - 1), 'Value', 1)
                        set(handles.slider8, 'Min', 1, 'Max', max(handles.fourth_dim_vector(:)), ... %change xlim
                            'SliderStep', [1, 1]/(max(handles.fourth_dim_vector(:)) - 1), 'Value', max(handles.fourth_dim_vector(:)))
                        guidata(hObject, handles)
                    case {'Cancel', ''}
                        return
                end
            end
        end
    end
end
hold(handles.axes1, 'off')
hold(handles.axes2, 'off')
guidata(hObject, handles)
%--------------------------------------------------------------------------------------------------------------------

%--------------------------------------------------------------------------------------------------------------------
% --- Executes on button press in SAVE DATA FROM AXES1.
function pushbutton9_Callback(hObject, eventdata, handles)
if ( strcmp(handles.save_format, 'Portable Network Graphics (*.png)') ||  strcmp(handles.save_format, 'JPEG Image (*.jpg)') )
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
    title(handles.external.(['ax' (num2str(handles.ext_plot))]), handles.axes1.Title.String, 'Interpreter', 'None')
else
    off_spring = findobj(handles.axes1, 'Type', 'Line');
    if ~isempty(off_spring)
        if (length(off_spring) > 1) && ( strcmp(handles.save_format, 'NIFTI (*.nii)') || strcmp(handles.save_format, 'DICOM (*.dcm)') )
            opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
            errordlg(['\fontsize{10}Can not save multiple Line objects to format ' handles.save_format], 'File save error', opts);
        else
            off_spring = off_spring.YData;
        end
    else
        off_spring = findobj(handles.axes1, 'Type', 'Image');
        off_spring = off_spring.CData;
    end
    switch handles.save_format
        case 'NIFTI (*.nii)'; file_filter = {'NIFTI *.nii'; 'DICOM *.dcm' ; 'MATLAB *.mat'; 'Text *.txt'; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat' };
        case 'DICOM (*.dcm)'; file_filter = {'DICOM *.dcm'; 'NIFTI *.nii' ; 'MATLAB *.mat'; 'Text *.txt'; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat'};
        case 'MATLAB (*.mat)'; file_filter = {'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii' ; 'Text *.txt' ; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat'};
        case 'Text (*.txt)'; file_filter = {'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii'; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat' };
        case 'DAT (*.dat)'; file_filter = {'DAT *.dat'; 'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii';  'Comma Separated Values *.csv'; 'Excel *.xlsx' };
        case 'Comma Separated Values (*.csv)'; file_filter = {'Comma Separated Values *.csv'; 'DAT *.dat'; 'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii'; 'Excel *.xlsx' };
        case 'Excel (*.xlsx)';  file_filter = {'Excel *.xlsx'; 'Comma Separated Values *.csv'; 'DAT *.dat'; 'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii' };
    end
    [file, path] = uiputfile(file_filter);
    if file ~= 0
        file_path = fullfile(path, file);
        try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
        save_data.this_folder = file_path; save_data.this_format = 0;
        save_data.this_name = 0; save_data.data_to_save = off_spring; save_this_file(save_data)
    end
end
guidata(hObject, handles)
%------------------------------------------------------------------------------------------------------------------

% --- Executes on button press in LOAD DATA ONTO AXES2.
function pushbutton10_Callback(hObject, eventdata, handles)
cla(handles.axes2, 'reset')
set(handles.axes2,'XTick', [], 'YTick', [])
set(handles.axes2, 'box', 'on')
handles.plotting_4d_data = false;
set(handles.slider1, 'Visible', 'on')
[file, path] = uigetfile('*');
if file == 0
else
    handles.data_onto_axes2 = true;
    path_to_axes2_data = fullfile(path, file);
    wrap_text = ['Visualising data from ' path_to_axes2_data];
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    [~, file_name, file_ext] = fileparts(path_to_axes2_data);
    %recognise files generated by DECONVOLVER and read them accordingly
    switch file_name
        case {'cbf_BzD', 'cbf_SVD', 'cbv_as_cbf_times_mtt_BzD', 'cbv_as_cbf_times_mtt_SVD', 'cbv_without_deconvolution', 'cmro2_BzD', 'cmro2_SVD',...
                'delay_BzD', 'delay_SVD', 'mtt_as_area_under_r_BzD', 'mtt_as_area_under_r_SVD', 'mtt_as_cbv_over_cbf_BzD', 'mtt_as_cbv_over_cbf_SVD',...
                'oef_BzD', 'oef_SVD', 'ttp_from_calculated_signal_SVD', 'ttp_from_fitted_signal_BzD', 'ttp_without_deconvolution'}; num_dims = 3;
        case {'residue_functions_BzD', 'residue_functions_SVD', 'control_points_BzD', 'disp_kernel_BzD'}; num_dims = 4;
        otherwise; num_dims = false;
    end
    axes2_pre_data = read_this_file(fullfile(path, file),num_dims);
    if isstruct(axes2_pre_data)
        axes2_data = squeeze(axes2_pre_data.the_data); %remove singular dimensions as they will confuse sliders
        handles.header_info = axes2_pre_data.the_info;
    else
        if isequal(axes2_pre_data, false)
            return
        else
            axes2_data = squeeze(axes2_pre_data);
        end
    end
    size_axes2_data = size(axes2_data);
    handles.size_axes2_data = size_axes2_data;
    handles.axes2_data = axes2_data;
    handles.path_to_axes2_data = path_to_axes2_data;
    handles.axes2_filename = file_name;
    clim_max = double(max(axes2_data(:))) + 1E-10;
    if clim_max <= 99; step_small = clim_max/1000; step_big = clim_max/999; end
    if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end
    if length(size_axes2_data) == 2 %2D or 1D
        if ~isempty(size_axes2_data(size_axes2_data == 1)) %1D data
            set(handles.slider2, 'Visible', 'off')
            set(handles.slider7, 'Visible', 'off')
            set(handles.slider11, 'Visible', 'off')
            set(handles.slider9, 'Visible', 'off')
            set(handles.pushbutton16, 'Visible', 'off')
            set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes2
            set(handles.togglebutton3, 'Visible', 'off')
            %                 wrap_text = 'Waiting for user input...';
            %                 set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
            opts.Interpreter = 'tex';
            prompt = {['\color{blue} \fontsize{10}You have loaded 1D data (',...
                num2str(size_axes2_data(1)), '\times', num2str(size_axes2_data(2)),') ',...
                ' onto axes1. Program will display  your data as a plot in the Cartesian plane.', ...
                sprintf(['\n' '\n']), '\color{blue} \fontsize{10}If applicable, provide x-axis sampling interval below.'...
                ' Leave edit field unchanged otherwise.']};
            defaults = {'1'};
            dims = [1 100];
            dlgtitle = 'Visualisation of 1D data';
            response = inputdlg(prompt, dlgtitle, dims, defaults, opts);
            if ~isempty(response) %Clicked OK
                d_x = str2double(response{1});
                x = 0:d_x:d_x*(max(size_axes2_data)-1);
                plot(x, axes2_data, 'k-','LineWidth',1, 'Parent', handles.axes2)
                title(handles.axes2, file_name, 'Interpreter', 'none')
            end
        else %if 2D image
            set(handles.slider2, 'Visible', 'off')
            set(handles.slider7, 'Visible', 'on')
            set(handles.slider11, 'Visible', 'on')
            set(handles.slider9, 'Visible', 'off')
            set(handles.pushbutton16, 'Visible', 'on')
            set(handles.pushbutton18, 'Visible', 'on')%rotate image on axes2
            set(handles.togglebutton3, 'Visible', 'on')
            handles.axes2_image = imshow(axes2_data, 'colormap', gray, 'Parent', handles.axes2);
            handles.axes2.CLim = [0 clim_max];
            title(handles.axes2, file_name, 'Interpreter', 'none')
            title(handles.axes2, file_name, 'Interpreter', 'none')
            set(handles.slider2, 'Visible', 'off') %slices
            set(handles.slider7, 'Visible', 'on') %color setting
            set(handles.slider7, 'Min', 10E-10, 'Max', clim_max, ...
                'SliderStep', [step_small, step_big], 'Value', clim_max)
            set(handles.slider11, 'Visible', 'on') %color setting
            set(handles.slider11, 'Min', 0, 'Max', clim_max-1E-10, ...
                'SliderStep', [step_small, step_big], 'Value', 0)
        end
    end
    
    if length(size_axes2_data) == 3 %3D data
        if ( strcmp(file_ext, '.png') || strcmp(file_ext, '.PNG') ...
                || strcmp(file_ext, '.jpg') || strcmp(file_ext, '.JPG') || ...
                strcmp(file_ext, '.JPEG') || strcmp(file_ext, '.gif') || strcmp(file_ext, '.ico'))
            handles.axes2_image = imshow(axes2_data, 'Parent', handles.axes2); %ensure correct visualisation of rgb photos
            set(handles.slider2, 'Visible', 'off')
            set(handles.slider7, 'Visible', 'on')
            set(handles.slider11, 'Visible', 'on')
            set(handles.slider9, 'Visible', 'off')
            set(handles.pushbutton16, 'Visible', 'on')
            set(handles.pushbutton18, 'Visible', 'on')%rotate image on axes2
            set(handles.togglebutton3, 'Visible', 'on')
            title(handles.axes2, file_name, 'Interpreter', 'none')
            guidata(hObject, handles)
        else
            set(handles.slider2, 'Visible', 'on')
            set(handles.slider7, 'Visible', 'on')
            set(handles.slider11, 'Visible', 'on')
            set(handles.slider9, 'Visible', 'off')
            set(handles.pushbutton16, 'Visible', 'on')
            set(handles.pushbutton18, 'Visible', 'on')%rotate image on axes2
            set(handles.togglebutton3, 'Visible', 'on')
            handles.axes2_image = imshow(squeeze(axes2_data(:,:,1)), 'colormap', gray, 'Parent', handles.axes2);
            handles.axes2.CLim = [0 clim_max];
            set(handles.slider2, 'Min', 1, 'Max', size_axes2_data(3), ... %slices
                'SliderStep', [1, 1]/(size_axes2_data(3) - 1), 'Value', 1)
            set(handles.slider7, 'Visible', 'on') %color setting
            set(handles.slider7, 'Min', 10E-10, 'Max', clim_max, ...
                'SliderStep', [step_small, step_big], 'Value', clim_max)
            set(handles.slider11, 'Visible', 'on') %color setting
            set(handles.slider11, 'Min', 0, 'Max', clim_max-1E-10, ...
                'SliderStep', [step_small, step_big], 'Value', 0)
            title(handles.axes2, [file_name, ': Slice 1'], 'Interpreter', 'none')
            guidata(hObject, handles)
        end
    end
    
    if length(size_axes2_data) == 4 %nasty
        %                 wrap_text = 'Waiting for user input...';
        %                 set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
        opts.Interpreter = 'tex';
        prompt = {['\color{blue} \fontsize{10}You have loaded 4D data (',...
            num2str(size_axes2_data(1)), '\times', num2str(size_axes2_data(2)), '\times', num2str(size_axes2_data(3)), ...
            '\times', num2str(size_axes2_data(4)), ') ',...
            'onto axes1. Program can either display your data as 1D plots in the Cartesian plane or image sequences.', ...
            ' Specifiy desired display (plots or images).'], [sprintf(['\n' '\n']), ...
            '\color{blue} \fontsize{10}If applicable, provide x-axis/temporal sampling interval below.'...
            ' Leave edit field unchanged otherwise.']};
        defaults = {'images', '1'};
        dims = [1 100];
        dlgtitle = 'Visualisation of 4D data';
        response = inputdlg(prompt, dlgtitle, dims, defaults, opts);
        if ~isempty(response) %Clicked OK
            d_t = str2double(response{2}); %extract temporal spacing; default is unity
            t = 0:d_t:d_t*(max(size_axes2_data(4))-1); %generate 4th dimension time vector
            handles.fourth_dim_vector = t;
            if ~strcmp(response{1}, 'plots') %if user selects image display mode
                %                     title(handles.axes1, file_name, 'Interpreter', 'none')
                set(handles.slider2, 'Visible', 'on')
                set(handles.slider7, 'Visible', 'on')
                set(handles.slider11, 'Visible', 'on')
                set(handles.slider9, 'Visible', 'on')
                set(handles.pushbutton16, 'Visible', 'on')
                set(handles.pushbutton18, 'Visible', 'on')%rotate image on axes2
                set(handles.togglebutton3, 'Visible', 'on')
                handles.axes2_image = imshow(squeeze(axes2_data(:,:,1,1)), 'colormap', gray, 'Parent', handles.axes2);
                handles.axes2.CLim = [0 clim_max];
                set(handles.slider2, 'Min', 1, 'Max', size_axes2_data(3), ... %1st adn 2nd dimensions
                    'SliderStep', [1, 1]/(size_axes2_data(3) - 1), 'Value', 1)
                dummy_data = axes2_data; dummy_data(~isfinite(dummy_data)) = []; dummy_data(isnan(dummy_data)) = [];
                clim_max = double(max(dummy_data(:))) + 1E-10; clear('dummy_data')
                if clim_max <= 99; step_small = clim_max/100; step_big = clim_max/99; end
                if clim_max > 99; step_small = 1/clim_max; step_big = step_small; end
                set(handles.slider7, 'Visible', 'on') %color setting
                set(handles.slider7, 'Min', 10E-10, 'Max', clim_max, ...
                    'SliderStep', [step_small, step_big], 'Value', clim_max)
                set(handles.slider11, 'Visible', 'on') %color setting
                set(handles.slider11, 'Min', 0, 'Max', clim_max-1E-10, ...
                    'SliderStep', [step_small, step_big], 'Value', 0)
                set(handles.slider9, 'Min', 1, 'Max', size_axes2_data(4), ... %4th dimension
                    'SliderStep', [1, 1]/(size_axes2_data(4) - 1), 'Value', 1)
                title(handles.axes2, [file_name, ': Slice 1 : Time-point 1'], 'Interpreter', 'none')
            end
            
            if strcmp(response{1}, 'plots') %if user selects plots display mode
                set(handles.slider2, 'Visible', 'off')
                set(handles.slider7, 'Visible', 'off')
                set(handles.slider11, 'Visible', 'off')
                set(handles.slider9, 'Visible', 'off')
                set(handles.pushbutton16, 'Visible', 'off')
                set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes2
                set(handles.togglebutton3, 'Visible', 'off')
                %                         wrap_text = 'Waiting for user input...';
                %                         set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
                opts.Interpreter = 'tex';
                num_plots = size_axes2_data(1)*size_axes2_data(2);
                prompt = {[sprintf('\n'), ...
                    '\color{blue} \fontsize{10}You are about to plot ',...
                    num2str(size_axes2_data(1)), '\times', num2str(size_axes2_data(2)), ' ( = ', num2str(num_plots), ...
                    ') plots on axes2. ',...
                    '\color{red} \fontsize{10}This is highly unrecommended.', ...
                    '\color{blue} \fontsize{10} Program allows selection of one or more ROIs from which to plot',...
                    '\color{blue} \fontsize{10} Press ''Define ROI'' to load image on which to place ROI(s).',...
                    sprintf(['\n' '\n']), ...
                    '\color{magenta} \fontsize{10}WARNING: Image will be loaded onto axes1. Axes1 will be cleared. ',...
                    sprintf('\n'), ...
                    '\color{blue} \fontsize{10}Select action.']}; %#ok<SPRINTFN>
                opts.Default = 'Define ROI';
                %                         dims = [1 100];
                dlgtitle = 'Visualisation of 4D data';
                response2 = questdlg(prompt, dlgtitle, 'Define ROI','Do not define ROI', 'Cancel', opts);
                switch response2 %user did not cancel dialog
                    case 'Define ROI' %user chose to define ROI
                        cla(handles.axes1, 'reset')
                        set(handles.axes1,'XTick', [], 'YTick', [])
                        set(handles.axes1, 'box', 'on')
                        roi_type = 'Freehand'; %type of ROI
                        [file, path] = uigetfile('*');
                        [~, file_name2, ~] = fileparts(fullfile(path, file));
                        %recognise files generated by DECONVOLVER and read them accordingly
                        switch file_name2
                            case {'cbf_BzD', 'cbf_SVD', 'cbv_as_cbf_times_mtt_BzD', 'cbv_as_cbf_times_mtt_SVD', 'cbv_without_deconvolution', 'cmro2_BzD', 'cmro2_SVD',...
                                    'control_points_BzD', 'delay_BzD', 'delay_SVD', 'mtt_as_area_under_r_BzD', 'mtt_as_area_under_r_SVD', 'mtt_as_cbv_over_cbf_BzD', 'mtt_as_cbv_over_cbf_SVD',...
                                    'oef_BzD', 'oef_SVD', 'ttp_from_calculated_signal_SVD', 'ttp_from_fitted_signal_BzD', 'ttp_without_deconvolution'}; num_dims = 3;
                            case {'residue_functions_BzD', 'residue_functions_SVD'}; num_dims = 4;
                            otherwise; num_dims = false;
                        end
                        if file == 0 %do nothing if no file is selected
                        else
                            axes1_pre_data = read_this_file(fullfile(path, file), num_dims);
                            if isstruct(axes1_pre_data)
                                axes1_data = squeeze(axes1_pre_data.the_data);
                                handles.header_info = axes1_pre_data.the_info;
                            else
                                if isequal(axes1_pre_data, false)
                                    return
                                else
                                    axes1_data = squeeze(axes1_pre_data);
                                end
                            end
                            size_axes1_data = size(axes1_data);
                            handles.size_axes1_data = size_axes1_data;
                            handles.axes1_data = axes1_data;
                            clim_max_2 = double(max(axes1_data(:))) + 1E-10;
                            if clim_max_2 <= 99; step_small = clim_max_2/100; step_big = clim_max_2/99; end
                            if clim_max_2 > 99; step_small = 1/clim_max_2; step_big = step_small; end
                            [~, file_name_2, file_ext_2] = fileparts(fullfile(path, file));
                            handles.axes1_filename = file_name_2;
                            handles.data_onto_axes1 = true;
                            guidata(hObject, handles)
                            if length(size_axes1_data) == 2 %loaded 2D image
                                set(handles.slider1, 'Visible', 'off')
                                set(handles.slider6, 'Visible', 'on')
                                set(handles.slider10, 'Visible', 'on')
                                set(handles.slider8, 'Visible', 'off')
                                set(handles.pushbutton15, 'Visible', 'on')
                                set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
                                set(handles.togglebutton2, 'Visible', 'on')
                                handles.axes2_roi_image = imshow(axes1_data, 'colormap', gray, 'Parent', handles.axes1);
                                handles.axes1.CLim = [0 clim_max_2];
                                set(handles.slider1, 'Visible', 'off') %slice change
                                set(handles.slider6, 'Visible', 'on') %color setting
                                set(handles.slider6, 'Min', 10E-10, 'Max', clim_max_2, ...
                                    'SliderStep', [step_small, step_big], 'Value', clim_max_2)
                                set(handles.slider10, 'Visible', 'on') %color setting
                                set(handles.slider10, 'Min', 0, 'Max', clim_max_2-1E-10, ...
                                    'SliderStep', [step_small, step_big], 'Value', 0)
                                title(handles.axes1, file_name_2, 'Interpreter', 'none')
                                guidata(hObject, handles)
                            end
                            
                            if length(size_axes1_data) == 3 %loaded 3D image
                                if ( strcmp(file_ext_2, '.png') || strcmp(file_ext_2, '.PNG') ...
                                        || strcmp(file_ext_2, '.jpg') || strcmp(file_ext_2, '.JPG') || ...
                                        strcmp(file_ext_2, '.JPEG') || strcmp(file_ext_2, '.gif') || strcmp(file_ext_2, '.ico'))
                                    handles.axes1_image = imshow(axes1_data, 'Parent', handles.axes1); %ensure correct visualisation of rgb photos
                                    set(handles.slider1, 'Visible', 'off')
                                    set(handles.slider6, 'Visible', 'on')
                                    set(handles.slider10, 'Visible', 'on')
                                    set(handles.slider8, 'Visible', 'off')
                                    set(handles.pushbutton15, 'Visible', 'on')
                                    set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
                                    set(handles.togglebutton2, 'Visible', 'on')
                                    title(handles.axes1, file_name_2, 'Interpreter', 'none')
                                    guidata(hObject, handles)
                                else
                                    set(handles.slider1, 'Visible', 'on')
                                    set(handles.slider6, 'Visible', 'on')
                                    set(handles.slider10, 'Visible', 'on')
                                    set(handles.slider8, 'Visible', 'off')
                                    set(handles.pushbutton15, 'Visible', 'on')
                                    set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
                                    set(handles.togglebutton2, 'Visible', 'on')
                                    handles.axes1_image = imshow(squeeze(axes1_data(:,:,1)), 'colormap', gray, 'Parent', handles.axes2);
                                    handles.axes1.CLim = [0 clim_max_2];
                                    set(handles.slider1, 'Min', 1, 'Max', size_axes1_data(3), ...
                                        'SliderStep', [1, 1]/(size_axes1_data(3) - 1), 'Value', 1)
                                    set(handles.slider6, 'Visible', 'on') %color setting
                                    set(handles.slider6, 'Min', 10E-10, 'Max', clim_max_2, ...
                                        'SliderStep', [step_small, step_big], 'Value', clim_max_2)
                                    set(handles.slider10, 'Visible', 'on') %color setting
                                    set(handles.slider10, 'Min', 0, 'Max', clim_max_2-1E-10, ...
                                        'SliderStep', [step_small, step_big], 'Value', 0)
                                    title(handles.axes1, [file_name_2, ': Slice 1'], 'Interpreter', 'none')
                                    guidata(hObject, handles)
                                end
                            end
                            
                            if length(size_axes1_data) == 4 %loaded 4D image
                                set(handles.slider1, 'Visible', 'on')
                                set(handles.slider6, 'Visible', 'on')
                                set(handles.slider10, 'Visible', 'on')
                                set(handles.slider8, 'Visible', 'on')
                                set(handles.pushbutton15, 'Visible', 'on')
                                set(handles.pushbutton17, 'Visible', 'off')%rotate image on axes1
                                set(handles.togglebutton2, 'Visible', 'on')
                                handles.axes1_image = imshow(squeeze(axes1_data(:,:,1,1)), 'colormap', gray, 'Parent', handles.axes1);
                                handles.axes1.CLim = [0 clim_max_2];
                                set(handles.slider1, 'Min', 1, 'Max', size_axes1_data(3), ... %1st adn 2nd dimensions
                                    'SliderStep', [1, 1]/(size_axes1_data(3) - 1), 'Value', 1)
                                dummy_data = axes1_data; dummy_data(~isfinite(dummy_data)) = []; dummy_data(isnan(dummy_data)) = [];
                                clim_max_2 = double(max(dummy_data(:))); clear('dummy_data')
                                if clim_max_2 <= 99; step_small = clim_max_2/100; step_big = clim_max_2/99; end
                                if clim_max_2 > 99; step_small = 1/clim_max_2; step_big = step_small; end
                                set(handles.slider6, 'Visible', 'on') %color setting
                                set(handles.slider6, 'Min', 10E-10, 'Max', clim_max_2, ...
                                    'SliderStep', [step_small, step_big], 'Value', clim_max_2)
                                set(handles.slider10, 'Visible', 'on') %color setting
                                set(handles.slider10, 'Min', 0, 'Max', clim_max_2-1E-10, ...
                                    'SliderStep', [step_small, step_big], 'Value', 0)
                                set(handles.slider8, 'Min', 1, 'Max', size_axes1_data(4), ... %4th dimension
                                    'SliderStep', [1, 1]/(size_axes1_data(4) - 1), 'Value', 1)
                                title(handles.axes1, [file_name_2, ': Slice 1 : Time-point 1'], 'Interpreter', 'none')
                                guidata(hObject, handles)
                            end
                            %OK, whatever image has been loaded is up
                            %and ready for ROI placement.
                            
                            %----------------------------------------------------------------------------------------------------------------
                            wrap_text = ['Click on the image to begin drawing of freehand ROI(s).  Double-click in ROI when done. '...
                                'Press any key to cancel.'...
                                'Use slider to change displayed slice.' ];
                            set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
                            blink_this(handles.edit7, 'b')
                            analyse_all = waitforbuttonpress;
                            n_rois = 128*128;
                            this_pos = 0;
                            handles.plot_colors = {'black', 'green', 'blue', 'red', 'cyan', 'magenta', 'yellow', rgb('Silver'), rgb('DarkGreen'), rgb('SkyBlue'), rgb('DarkRed'), rgb('CadetBlue'),...
                                rgb('Indigo'), rgb('Brown'), rgb('Gray'), rgb('HotPink'), rgb('DarkOrange'), rgb('DodgerBlue'), rgb('RosyBrown'), rgb('Lime'), rgb('OrangeRed'), rgb('GreenYellow'),...
                                rgb('Purple'), rgb('Chocolate'), rgb('SaddleBrown'), rgb('DarkKhaki'), rgb('Pink'), rgb('Salmon'), rgb('Olive'), rgb('LightSteelBlue')};
                            if ~analyse_all %mouse clicked
                                for this_roi = 1:n_rois
                                    if this_roi > 1
                                        wrap_text = ['Click to add another ROI. Press: Shift + [R (rectangle), C (circle), E (ellipse)].',...
                                            'Press Esc to end.', ' Press Spacebar to reset axes.'];
                                        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
                                        enough = waitforbuttonpress;
                                        this_key = double(get(handles.figure1, 'CurrentCharacter'));
                                        if enough %if user presses key
                                            switch this_key
                                                case 27; break %Esc
                                                case 82; roi_type = 'Rectangle'; %Shift R
                                                    %                    case 80; roi_type = 'Point'; %Shift P
                                                case 67; roi_type = 'Circle'; %Shift C
                                                case 69; roi_type = 'Ellipse'; %Shift E
                                                case 70; roi_type = 'Freehand'; %Shift F
                                                case 32 %Spacebar
                                                    current_title = handles.axes1.Title.String;
                                                    cla(handles.axes2, 'reset')
                                                    set(handles.axes2,'XTick', [], 'YTick', [])
                                                    set(handles.axes2, 'box', 'on')
                                                    handles.axes2.Title.String = current_title;
                                                    for rr = 1:this_pos
                                                        if ishandle(handles.fourd_rois.(['roi_' num2str(rr)]))
                                                            delete(handles.fourd_rois.(['roi_' num2str(rr)]))
                                                            curves = fieldnames(handles.(['roi_curves' num2str(rr)]));
                                                            for cc = 1:numel(curves)
                                                                delete(handles.(['roi_curves' num2str(rr)]).(['curve_' num2str(cc)]))
                                                            end
                                                        end
                                                    end
                                                    this_pos = 0;
                                            end
                                        end
                                    end
                                    
                                    wrap_text = ['Select region to analyse. Click on the image to begin drawing of ROI.  Double-click in ROI when done. '...
                                        'Use slider to change displayed slice.' ];
                                    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
                                    switch roi_type
                                        case 'Freehand'; fourd_roi = drawfreehand(handles.axes1, 'Color', 'yellow', 'LineWidth', 0.1);
                                        case 'Rectangle'; fourd_roi = drawrectangle(handles.axes1, 'Color', 'yellow', 'LineWidth', 0.1);
                                        case 'Circle'; fourd_roi = drawcircle(handles.axes1, 'Color', 'yellow', 'LineWidth', 0.1);
                                        case 'Ellipse'; fourd_roi = drawellipse(handles.axes1, 'Color', 'yellow', 'LineWidth', 0.1);
                                    end
                                    fourd_roi_pos = customWait(fourd_roi);
                                    this_pos = this_pos + 1;
                                    if ishandle(fourd_roi); set(fourd_roi, 'Color', handles.plot_colors{this_pos}); end
                                    handles.fourd_rois.(['roi_' num2str(this_pos)]) = fourd_roi;
                                    drawn_slice = round(get(handles.slider1, 'Value'));
                                    
                                    %if user chose not to generate mask, create a mask true everywhere
                                    curve = 0;
                                    for x = 1:handles.size_axes2_data(1)
                                        for y = 1:handles.size_axes2_data(2) %for each pixel
                                            if inROI(handles.fourd_rois.(['roi_' num2str(this_pos)]), x, y) %if pixel exists in at least one of those rois
                                                curve = curve + 1;
                                                handles.(['roi_curves' num2str(this_pos)]).(['curve_' num2str(curve)]) = plot(t, squeeze(axes2_data(y,x,drawn_slice,:)),...
                                                    'Color', handles.plot_colors{this_pos}, 'Parent', handles.axes2);
                                                hold(handles.axes2, 'on')
                                            end
                                        end
                                    end
                                    title(handles.axes2, file_name, 'Interpreter', 'none')
                                end
                            else %if key pressed, do nothing
                            end
                        end
                        %----------------------------------------------------------------------------------------------------------------
                    case 'Do not define ROI' %user edited the roi selection field = does not want to select ROI
                        handles.this_object = hObject;
                        set(handles.slider2, 'Visible', 'on')
                        set(handles.slider7, 'Visible', 'off')
                        set(handles.slider11, 'Visible', 'off')
                        set(handles.slider9, 'Visible', 'off')
                        set(handles.pushbutton16, 'Visible', 'off')
                        set(handles.pushbutton18, 'Visible', 'off')%rotate image on axes2
                        set(handles.togglebutton3, 'Visible', 'off')
                        title(handles.axes2, [file_name ': Slice 1'], 'Interpreter', 'none')
                        pause(0.5)
                        plot_4d_data(handles, axes2_data, 1, handles.axes2, file_name)
                        handles.plotting_4d_data = true; %Here
                        set(handles.slider2, 'Min', 1, 'Max', size_axes2_data(3), ... %1st adn 2nd dimensions
                            'SliderStep', [1, 1]/(size_axes2_data(3) - 1), 'Value', 1)
                        set(handles.slider9, 'Min', 1, 'Max', max(handles.fourth_dim_vector(:)), ... %change xlim
                            'SliderStep', [1, 1]/(max(handles.fourth_dim_vector(:)) - 1), 'Value', max(handles.fourth_dim_vector(:)))
                        guidata(hObject, handles)
                    case {'Cancel', ''}
                        return
                end
            end
        end
    end
end
hold(handles.axes2, 'off')
hold(handles.axes1, 'off')
guidata(hObject, handles)
%---------------------------------------------------------------------------------------------------------

% --- Executes on button press in SAVE DATA FROM AXES2.
function pushbutton11_Callback(hObject, eventdata, handles)
if ( strcmp(handles.save_format, 'Portable Network Graphics (*.png)') ||  strcmp(handles.save_format, 'JPEG Image (*.jpg)') )
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
else
    off_spring = findobj(handles.axes2, 'Type', 'Line');%handles.axes1.Children
    if ~isempty(off_spring)
        if (length(off_spring) > 1) && ( strcmp(handles.save_format, 'NIFTI (*.nii)') || strcmp(handles.save_format, 'DICOM (*.dcm)') )
            wrap_text = ['Cannot save multiple Line objects to format: ' handles.save_format];
            set(handles.edit7, 'String', wrap_text, 'ForeGroundColor', 'r')
        else
            off_spring = off_spring.YData;
        end
    else
        off_spring = findobj(handles.axes2, 'Type', 'Image');
        off_spring = off_spring.CData;
    end
    switch handles.save_format
        case 'NIFTI (*.nii)'; file_filter = {'NIFTI *.nii'; 'DICOM *.dcm' ; 'MATLAB *.mat'; 'Text *.txt'; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat' };
        case 'DICOM (*.dcm)'; file_filter = {'DICOM *.dcm'; 'NIFTI *.nii' ; 'MATLAB *.mat'; 'Text *.txt'; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat'};
        case 'MATLAB (*.mat)'; file_filter = {'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii' ; 'Text *.txt' ; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat'};
        case 'Text (*.txt)'; file_filter = {'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii'; 'Comma Separated Values *.csv'; 'Excel *.xlsx'; 'DAT *.dat' };
        case 'DAT (*.dat)'; file_filter = {'DAT *.dat'; 'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii';  'Comma Separated Values *.csv'; 'Excel *.xlsx' };
        case 'Comma Separated Values (*.csv)'; file_filter = {'Comma Separated Values *.csv'; 'DAT *.dat'; 'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii'; 'Excel *.xlsx' };
        case 'Excel (*.xlsx)';  file_filter = {'Excel *.xlsx'; 'Comma Separated Values *.csv'; 'DAT *.dat'; 'Text *.txt'; 'MATLAB *.mat'; 'DICOM *.dcm'; 'NIFTI *.nii' };
    end
    [file, path] = uiputfile(file_filter);
    if file ~= 0
        file_path = fullfile(path, file);
        try save_data.header_info = handles.header_info; catch; save_data.header_info = []; end
        save_data.this_folder = file_path; save_data.this_format = 0;
        save_data.this_name = 0; save_data.data_to_save = off_spring; save_this_file(save_data)
    end
end
guidata(hObject, handles)

% --- Executes on button press in ADDITIONAL CALCULATIONS: OEF.
function checkbox26_Callback(hObject, eventdata, handles)
handles.calculate_oef = get(hObject,'Value');
guidata(hObject, handles)


% --- Executes on button press in ADDITIONAL CALCULATIONS: CMRO2.
function checkbox27_Callback(hObject, eventdata, handles)
handles.calculate_cmro2 = get(hObject,'Value');
guidata(hObject, handles)


% --- Executes on button press in ADDITIONAL CALCULATIONS: R10.
function checkbox28_Callback(hObject, eventdata, handles)
handles.calculate_r10 = get(hObject,'Value');
guidata(hObject, handles)


% --- Executes on button press in ADDITIONAL CALCULATIONS: R50..
function checkbox29_Callback(hObject, eventdata, handles)
handles.calculate_r50 = get(hObject,'Value');
guidata(hObject, handles)


% --- Executes on selection change in DATA DISPLAY OPTIONS FOR AXES1.
function popupmenu11_Callback(hObject, eventdata, handles)
display_choices_axes2 = cellstr(get(handles.popupmenu12,'String'));
display_choice_axes2 = display_choices_axes2{get(handles.popupmenu12,'Value')};
handles.show_cbf = false;handles.show_mtt = false;handles.show_oef = false;handles.show_cmro2 = false;
handles.show_r10 = false; handles.show_r50 = false; handles.plot_residue_funcs = false; handles.show_delay = false;
switch display_choice_axes2
    case 'CBF'; handles.show_cbf = true;
    case 'MTT'; handles.show_mtt = true;
    case 'OEF'; handles.show_oef = true;
    case 'CMRO2'; handles.show_cmro2 = true;
    case 'R10'; handles.show_r10 = true;
    case 'R50';  handles.show_r50 = true;
    case 'Residue functions'; handles.plot_residue_funcs = true;
    case 'Delay'; handles.show_delay = true;
    case 'Nothing'; handles.axes2_vaccant = true;
end
display_choices = cellstr(get(hObject,'String'));
display_choice = display_choices{get(hObject,'Value')};

switch display_choice
    case 'CBF'; handles.show_cbf = true; handles.cbf_axes = handles.axes1; handles.cbf_slider = handles.slider1; handles.cbf_window_slider_up = handles.slider6; handles.cbf_window_slider_low = handles.slider10;
    case 'MTT'; handles.show_mtt = true; handles.mtt_axes = handles.axes1; handles.mtt_slider = handles.slider1; handles.mtt_window_slider_up = handles.slider6; handles.mtt_window_slider_low = handles.slider10;
    case 'OEF'; handles.show_oef = true; handles.oef_axes = handles.axes1; handles.oef_slider = handles.slider1; handles.oef_window_slider_up = handles.slider6; handles.oef_window_slider_low = handles.slider10;
    case 'CMRO2'; handles.show_cmro2 = true; handles.cmro2_axes = handles.axes1; handles.cmro2_slider = handles.slider1; handles.cmro2_window_slider_up = handles.slider6; handles.cmro2_window_slider_low = handles.slider10;
    case 'R10'; handles.show_r10 = true; handles.r10_axes = handles.axes1; handles.r10_slider = handles.slider1; handles.r10_window_slider_up = handles.slider6; handles.r10_window_slider_low = handles.slider10;
    case 'R50';  handles.show_r50 = true; handles.r50_axes = handles.axes1; handles.r50_slider = handles.slider1; handles.r50_window_slider_up = handles.slider6; handles.r50_window_slider_low = handles.slider10;
    case 'Residue functions'; handles.plot_residue_funcs = true; handles.rt_axes = handles.axes1; handles.rt_slider = handles.slider1; handles.rt_window_slider_up = handles.slider6; handles.rt_window_slider_low = handles.slider10;
    case 'Delay'; handles.show_delay = true; handles.delay_axes = handles.axes1; handles.delay_slider = handles.slider1; handles.delay_window_slider_up = handles.slider6; handles.delay_window_slider_low = handles.slider10;
    case 'Nothing'; handles.axes1_vaccant = true;
end
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in  DATA DISPLAY OPTIONS FOR AXES2
function popupmenu12_Callback(hObject, eventdata, handles)
display_choices_axes1 = cellstr(get(handles.popupmenu11,'String'));
display_choice_axes1 = display_choices_axes1{get(handles.popupmenu11,'Value')};
handles.show_cbf = false;handles.show_mtt = false;handles.show_oef = false;handles.show_cmro2 = false;
handles.show_r10 = false; handles.show_r50 = false;handles.plot_residue_funcs = false; handles.show_delay = false;
switch display_choice_axes1
    case 'CBF'; handles.show_cbf = true;
    case 'MTT'; handles.show_mtt = true;
    case 'OEF'; handles.show_oef = true;
    case 'CMRO2'; handles.show_cmro2 = true;
    case 'R10'; handles.show_r10 = true;
    case 'R50';  handles.show_r50 = true;
    case 'Delay'; handles.show_delay = true;
    case 'Residue functions'; handles.plot_residue_funcs = true;
    case 'Nothing'; handles.axes1_vaccant = true;
end
display_choices = cellstr(get(hObject,'String'));
display_choice = display_choices{get(hObject,'Value')};
switch display_choice
    case 'CBF'; handles.show_cbf = true; handles.cbf_axes = handles.axes2; handles.cbf_slider = handles.slider2; handles.cbf_window_slider_up = handles.slider7;  handles.cbf_window_slider_low = handles.slider11;
    case 'MTT'; handles.show_mtt = true; handles.mtt_axes = handles.axes2; handles.mtt_slider = handles.slider2;  handles.mtt_window_slider_up = handles.slider7; handles.mtt_window_slider_low = handles.slider11;
    case 'OEF'; handles.show_oef = true; handles.oef_axes = handles.axes2; handles.oef_slider = handles.slider2;  handles.oef_window_slider_up = handles.slider7;  handles.oef_window_slider_low = handles.slider11;
    case 'CMRO2'; handles.show_cmro2 = true; handles.cmro2_axes = handles.axes2; handles.cmro2_slider = handles.slider2;  handles.cmro2_window_slider_up = handles.slider7; handles.cmro2_window_slider_low = handles.slider11;
    case 'R10'; handles.show_r10 = true; handles.r10_axes = handles.axes2; handles.r10_slider = handles.slider2;  handles.r10_window_slider_up = handles.slider7; handles.r10_window_slider_low = handles.slider11;
    case 'R50';  handles.show_r50 = true; handles.r50_axes = handles.axes2; handles.r50_slider = handles.slider2;  handles.r50_window_slider_up = handles.slider7;  handles.r50_window_slider_low = handles.slider11;
    case 'Residue functions'; handles.plot_residue_funcs = true; handles.rt_axes = handles.axes2; handles.rt_slider = handles.slider2; handles.rt_window_slider_up = handles.slider7; handles.rt_window_slider_low = handles.slider11;
    case 'Delay'; handles.show_delay = true; handles.delay_axes = handles.axes2; handles.delay_slider = handles.slider2; handles.delay_window_slider_up = handles.slider7; handles.delay_window_slider_low = handles.slider11;
    case 'Nothing'; handles.axes2_vaccant = true;
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function popupmenu12_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in CONTINUE WITH SELECTED SLICE
function pushbutton14_Callback(hObject, eventdata, handles)
uiresume(handles.figure1)
guidata(hObject, handles)


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
function slider8_Callback(hObject, eventdata, handles)
time_point = round(get(hObject,'Value'));
current_slice = round(get(handles.slider1,'Value'));
if handles.plotting_4d_data
    xlim(handles.axes1, [0 time_point])
elseif handles.data_is_for_aif
    set(handles.aif_conc_image, 'CData', squeeze(handles.aif_conc_data(:,:,current_slice, time_point))  )
    title(handles.axes1, ['Concentration' ': Slice ' num2str(current_slice) ' : Time-point ' num2str(time_point)],...
        'Interpreter', 'None')
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
