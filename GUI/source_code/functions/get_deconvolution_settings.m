function handles_update = get_deconvolution_settings(handles, caller)
%        This function accepts the handles structure and the variable
%        'caller' (which shows which deconvolution algorithm is in use).
%        The function prompts the user for relevant settings (the
%        hematocrit kappa, OI, Psvd and range of deconvolution (slice
%        range)). An updated version of the handles structure is the
%        function's output.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     wrap_text = 'Waiting for user input...'; %display this message
     set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
     try %try to extract already-provided slice range
        current_z1 = num2str(handles.slice_range(5)); current_z2 = num2str(handles.slice_range(6));
     catch %otherwise use all slices as default slice range
        current_z1 = num2str(1); current_z2 = num2str(handles.img_size(3));
     end
     %try to get already provided kappa, Psvd and OI, otherwise use the
     %default values below.
    try current_kappa = num2str(handles.kappa); catch; current_kappa = '0.705'; end   
    try current_psvd = num2str(handles.Psvd); catch; current_psvd = '0.2'; end 
    try current_OIndex = num2str(handles.OIndex); catch; current_OIndex = '0.095'; end
    %______________________________________________________________________
    defaults = {current_kappa, current_psvd, current_OIndex, current_z1, current_z2};  %default settings    
    opts.Interpreter = 'tex'; %use the tex interpreter to change font size and colour
    prompt = {'\color{blue} \fontsize{10}Hematocrit constant ({\kappa_H}) [cm{^3}/g] :', ...
    '\color{blue} \fontsize{10}Psvd:', ...
    '\color{blue} \fontsize{10}Oscillation index (OI):', ...
        ['\color{blue} \fontsize{10}Deconvolution slice range:                                           '...
        '\color{blue} \fontsize{10}From slice number:'],...
        '\color{blue} \fontsize{10}To slice number:'}; %prompt to be displayed on input dialog
    switch caller %decide which settings to ask for depending on which deconvolution algorithm is in use (caller)
        case 'bzd'
            prompt(2:3) = [];
            defaults(2:3) = [];
            dlgtitle = 'DECONVOLVER: Edit settings for BzD';
            sl_indices = 2:3; OIndex_index = false; psvd_index = false;
        case 'osvd'
            prompt(2) = [];
            defaults(2) = [];
            dlgtitle = 'DECONVOLVER: Edit settings for oSVD';
            sl_indices = 3:4;
            OIndex_index = 2; psvd_index = false;
        case 'csvd'
            prompt(3) = [];
            defaults(3) = [];
            dlgtitle = 'DECONVOLVER: Edit settings for cSVD';
            sl_indices = 3:4; OIndex_index = false; psvd_index = 2;
        case 'ssvd'
            prompt(3) = [];
            defaults(3) = [];
            dlgtitle = 'DECONVOLVER: Edit settings for sSVD';
            sl_indices = 3:4; OIndex_index = false; psvd_index = 2;
    end
    %______________________________________________________________________
    dims = [1 80]; %dimensions of input dialog window
    response = inputdlg(prompt, dlgtitle, dims, defaults, opts); %invoke input dialog
    %______________________________________________________________________
    if isempty(response) %if use clicks cancel, report and terminate
        wrap_text = 'Analysis terminated by user. Re-run to select data.';
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
        handles.user_cancelled = true; %let main program know that user cancelled the dialog box
        handles_update = handles;
        return
    end
    %______________________________________________________________________
    handles.kappa = str2double(response{1}); %otherwise get kappa
    if OIndex_index; handles.OIndex = str2double(response{OIndex_index}); end %get OI
    if psvd_index; handles.Psvd = str2double(response{psvd_index}); end %get Psvd
    z1 = str2num(response{sl_indices(1)});  z2 = str2num(response{sl_indices(2)}); %get slice range [z1 z2]
    handles.slice_range = [1, handles.img_size(1), 1, handles.img_size(2), z1, z2];%deconvolution functions 'want' the slice range in the form [1 128 1 128 z1 z2]
    handles.user_cancelled = false; %user did not cancel
    handles_update = handles; %update handles structure
   end