function [te, tr] = get_te_tr(handles)
%       This function accepts the handles structure and returns Echo time
%       and Repetition time, either taken from the header of the input data
%       or requested from the user.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	path_to_data = handles.path_to_data;
	if isfile(path_to_data) %means user chose nifti or dicom 4d volume
        [~, ~, data_fmt] = fileparts(path_to_data);
    else
        if isfolder(path_to_data) %if path is a folder
            folder_contents = dir(path_to_data);
            file_names = {folder_contents.name};
            file_to_use = fullfile(path_to_data, file_names{round(length(file_names)/2)});
            path_to_data = file_to_use;
            data_fmt = '.dcm';
        end
    end
    
    switch data_fmt
		case {'.nii', '.gz'}
            opts.Interpreter = 'tex';
            opts.Default = 'Show header';
            response = questdlg(['\fontsize{10}Program cannot extract TE/TR from the selected data. '...
                'Program will attempt to display the header file for the selected data and request manual input of these parameters by the user. Select action.'], ...
            'Reading NIFTI file','Show header', 'Do not show header', 'Terminate', opts);
            switch response
                case 'Show header'
                    try delete('header.txt'); catch; end
                    HEADER = niftiinfo(path_to_data);
                    diary header.txt
                    disp(HEADER)
                    diary off
                    try 
                        winopen('header.txt')
                    catch
                        try
                            open('header.txt')
                        catch
                            opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
                            errordlg(['\fontsize{10}Failed to either read or open header file.'], 'Reading header file', opts);
                        end
                    end
                case 'Do not show header'
                case 'Terminate'; te = false; tr = false; return
                case ''; te = false; tr = false; return
            end
        case '.dcm'
            try
                HEADER = dicominfo(path_to_data);
                te = HEADER.EchoTime/1E3;
                tr = HEADER.RepetitionTime/1E3;
                handles.te = te;
                handles.tr = tr;
                opts.Interpreter = 'tex';
                opts.Default = 'Show header';
                response = questdlg(['\fontsize{10}TE/TR were successfully extracted from the selected data. '...
                'Program will attempt to display the header file for the selected data and request verification of these parameters by the user. Select action.'], ...
                'Reading DICOM file','Show header', 'Do not show header', 'Terminate', opts);
                switch response
                    case 'Show header'
                        try delete('header.txt'); catch; end
                        diary header.txt
                        disp(HEADER)
                        diary off
                        try 
                            winopen('header.txt')
                        catch
                            try
                                open('header.txt')
                            catch
                                opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
                                errordlg(['\fontsize{10}Failed to either read or open header file.'], 'Reading header file', opts);
                            end
                        end
                    case    'Do not show header'
                    case 'Terminate'; te = false; tr = false; return
                    case ''; te = false; tr = false;   return
                end
            catch %did not succeed in reading te tr
                opts.Interpreter = 'tex';
                opts.Default = 'Show header';
                response = questdlg(['\fontsize{10}Program could not extract from the selected data. '...
                'Program will attempt to display the header file for the selected data and request verification of these parameters by the user. Select action.'], ...
                'Reading DICOM file','Show header', 'Do not show header', 'Terminate', opts);
                switch response
                    case 'Show header'
                        try delete('header.txt'); catch; end
                        HEADER = dicominfo(path_to_data);
                        diary header.txt
                        disp(HEADER)
                        diary off
                        try 
                            winopen('header.txt')
                        catch
                            try
                                open('header.txt')
                            catch
                                opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
                                errordlg(['\fontsize{10}Failed to either read or open header file.'], 'Reading header file', opts);
                            end
                        end
                    case    'Do not show header'
                    case 'Terminate'; te = false; tr = false; return
                    case ''; te = false; tr = false;   return
                end
            end
     end

    wrap_text = 'Waiting for user input...';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
    opts.Interpreter = 'tex';
    prompt = {'\color{blue} \fontsize{10} TE [ms] :', '\color{blue} \fontsize{10} TR [ms] :'};
    dlgtitle = 'DECONVOLVER: Edit input data setup';
    dims = [1 80];
    try
        current_te = num2str(handles.te*1E3); 
        current_tr = num2str(handles.tr*1E3); 
    catch
        current_te = '';
        current_tr = '';
    end
    
    defaults = {current_te, current_tr};
    
    response = inputdlg(prompt, dlgtitle, dims, defaults, opts);
    
    if isempty(response) %if use clicks cancel, report and terminate
        wrap_text = 'Analysis terminated by user. Re-run to select data.';
        set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
        te = false; tr = false;
        return
    end 
    
    te = str2double(response{1})*1E-3;
    tr = str2double(response{2})*1E-3;
end
