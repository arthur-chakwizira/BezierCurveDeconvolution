function out_data = read_this_file(file_path, expected_dims)
%        This function accepts a path to a file and reads the data in the
%        file. expected_dims is the expected number of dimensions. If false,
%        the function makes no restrictions on the output data size. The
%        output of the function (out_data) is the data in the input file.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        
        [path_to_folder, file_name, file_extension] = fileparts(file_path);
        
        if expected_dims == false
            restricted_size = false; %do not restrict output data size if expected dimensions are not specified
        else
            restricted_size = true;
        end
        
        format_match = false; %this will remain false if the loaded data does not match any of the recognised file formats
        
        %NIFTI
        %------------------------------------------------------------------
        if isfile(file_path) && ( strcmp(file_extension, '.nii') || strcmp(file_extension, '.nii.gz') || strcmp(file_extension, '.gz'))
            format_match = true; %update format_match if format matches nifti
            out_data = niftiread(file_path); %read data from file
        end
        %------------------------------------------------------------------        
        
        %DICOM (IMA speculative right now; will be implemented later)
        %------------------------------------------------------------------
        %Single file (could be multiframe dicom or just one dicom image)
        if isfile(file_path) && (strcmp(file_extension, '.dcm') ||  strcmp(file_extension, '.ima')) %if dicom file
            format_match = true; %update format match
            dicom_data = dicomread(file_path); %read dicom file
            dicom_info = get_dicom_header(file_path); %extract header info
            out_data.the_data = dicom_data; %save volume to struct out_data
            out_data.the_info = dicom_info; %save header to struct out_data
            if restricted_size %if expected dimensions specified 
                if ndims(out_data) ~= expected_dims %if read data does not match expectations
%                     try
                        folder_contents = dir(path_to_folder); %get paths in folder
                        file_names = {folder_contents.name}; %get file names of Dicom series
                        fl = length(file_names);
                        to_be_removed = NaN(fl,1); %indices of invalid files
                        remove_idx = 0;
                        p_bar = waitbar(0, ['Verifying file validity: file 0 of ' , num2str(fl)], 'Name', 'Reading Dicom series');
                        for f = 1:fl
                            if isgraphics(p_bar)
                                waitbar(f/fl, p_bar, ['Verifying file validity: file ' num2str(f) ' of ' num2str(fl)]) %update progress bar
                            else
                                out_data = false; %return if user closes progress bar
                                return
                            end
                            [~, ~, exten] = fileparts(file_names{f}); %extract file extension for every file
                            if ~strcmp(exten, '.dcm'); remove_idx = remove_idx + 1; to_be_removed(remove_idx) = f; end %mark it for removal if it is not a dicom file
                        end
                        to_be_removed(isnan(to_be_removed)) = [];
                        file_names(to_be_removed) = []; %remove invalid files
                        first_file_info = dicominfo(fullfile(path_to_folder, file_names{1})); %get header info of first file in series
                        first_slice_locn = first_file_info.SliceLocation; %get the corresponding slice location
                        slice_spacing = first_file_info.SpacingBetweenSlices; %get slice spacing
                        last_file_info = dicominfo(fullfile(path_to_folder, file_names{end})); %get header info for last file in series
                        last_slice_locn = last_file_info.SliceLocation; %and corresponding slice location
                        num_slices = round(abs(first_slice_locn-last_slice_locn)/slice_spacing + 1); %determine number of slices
                        num_time_points = first_file_info.NumberOfTemporalPositions; %and number of time points
                        num_rows = double(first_file_info.Rows); %number of rows
                        num_columns = double(first_file_info.Columns); %columnes
                        all_slice_locs = NaN(num_slices, 1); 
                        new_slice = 0; %will change value when slice changes
                        if ishandle(p_bar); delete(p_bar); end
                        len_fn = length(file_names);
                        p_bar = waitbar(0, ['Preparing: file 0 of ' , num2str(len_fn)], 'Name', 'Reading Dicom series');
                        for ff = 1:len_fn
                            if isgraphics(p_bar)
                                waitbar(ff/len_fn, p_bar, ['Preparing: file ' num2str(ff) ' of ' num2str(len_fn)])
                            else
                                out_data = false;
                                return
                            end
                            tmp_tmp_info = dicominfo(fullfile(path_to_folder, file_names{ff})); %get info of every file
                            all_file_info.(['file_' num2str(ff)]) = tmp_tmp_info; %save info
                            tmp_slc_loc = tmp_tmp_info.SliceLocation; %get slice location
                            if ~ismember(tmp_slc_loc, all_slice_locs) %if slice has not been encountered before, increase slice counter
                                new_slice = new_slice+1;
                                all_slice_locs(new_slice) = tmp_slc_loc; %gather all slice locations
                            end
                        end
                        if ishandle(p_bar); delete(p_bar); end
                        p_bar = waitbar(0, ['Reading: file 0 of ' , num2str(len_fn)], 'Name', 'Reading Dicom series');                        
                        dicom_data = zeros(num_rows, num_columns, num_slices, num_time_points);
                                for fil = 1:len_fn
                                    if isgraphics(p_bar)
                                        waitbar(fil/length(file_names), p_bar, ['Reading: file ' num2str(fil) ' of ' num2str(len_fn)])
                                    else
                                        out_data = false;
                                        return
                                    end                                    
                                    tmp_data = dicomread(fullfile(path_to_folder, file_names{fil}));
                                    tmp_info = all_file_info.(['file_' num2str(fil)]);
                                    for slc = 1:num_slices
                                        slc_loc = all_slice_locs(slc);
                                        if tmp_info.SliceLocation == slc_loc
                                            for tim = 1:num_time_points
                                                if tmp_info.TemporalPositionIdentifier == tim
                                                    dicom_data(:,:,slc,tim) = tmp_data;
                                                end
                                            end
                                        end
                                    end
                                end
                        if ishandle(p_bar); delete(p_bar); end
%                         if ishandle(mbox); delete(mbox); end %delete the message box if user has not closed it yet
                        out_data.the_data = dicom_data; %save dicom data and
                        out_data.the_info = dicom_info; %header info
                        if ndims(out_data.the_data) == expected_dims %if the data matches required dimensions
                            return %we are good
                        end
%                     catch
%                         out_data = false; %if reading from folder does not work, report failur
%                     end
                else %if the data in the first place matches specified dimensions, we are done
                    return
                end
            end
        end
        
        %Volume/folder
        if isfolder(file_path) %user selected a folder
            try
                [dicom_data, ~, ~] = dicomreadVolume(path_to_folder); %try to read the volume
                dicom_info = get_dicom_header(path_to_folder); %get header info
                out_data.the_data = dicom_data; %save data
                out_data.the_info = dicom_info; %save info
                format_match = true;
            catch
                out_data = false;
            end
        end
        %------------------------------------------------------------------
        
        
        %PNG, JPEG, GIF
        %------------------------------------------------------------------
        if isfile(file_path) && ( strcmp(file_extension, '.png') || strcmp(file_extension, '.PNG') ...
                || strcmp(file_extension, '.jpg') || strcmp(file_extension, '.JPG') || ...
                strcmp(file_extension, '.JPEG') || strcmp(file_extension, '.gif') || strcmp(file_extension, '.ico'))
            format_match = true; %update format match if file extension matches any of the above
            out_data = imread(file_path); %read the file
        end
        %------------------------------------------------------------------
        
        %EXCEL: xlsx, xls, xlsm, xltx, xltm,...
        %------------------------------------------------------------------
        if isfile(file_path) && ( strcmp(file_extension, '.xlsx') || strcmp(file_extension, '.xls') || strcmp(file_extension, '.xlsm') || strcmp(file_extension, '.xltx') || strcmp(file_extension, '.xltm'))
            format_match = true; %update format match if file is Excel file
            [~, sheets, ~ ] = xlsfinfo(file_path); %get sheet names
            if length(sheets) > 1 %if there is more than 1 sheet
                %must ask user from which sheet to read
                %__________________________________________________________
               opts.Interpreter = 'tex';
                opts.Default = 'Do not start Excel';
                response = questdlg(['\fontsize{10}The selected Spreadsheet contains ' num2str(length(sheets)) ' sheets. Program will request sheet name to read data from. Select action.'], ...
                'Reading Excel file','Start Excel', 'Do not start Excel', 'Cancel', opts); %invoke this question dialog
                switch response
                        case 'Start Excel'; open_excel = true; %user pressed the Start Excel button
                        case 'Do not start Excel'; open_excel = false; %user pressed the Do not start Excel button
                        case 'Cancel'; out_data = false; return %user pressed cancel
                        case ''; out_data = false;  return %user closed the question dialog
                end
                if open_excel; try winopen(file_path); pause(1); catch; open(file_path); pause(1); end; end %try open Excel file in Excel, otherwise in MATLAB
                prompt = {'\color{blue} \fontsize{10}Enter sheet name', '\color{blue} \fontsize{10}Enter cell range ( example A1:A40 ). Leave edit field unchanged to read entire sheet.'};
                dlgtitle = 'Specify data location in Excel file';  %invoke input dialog to ask user for sheet name and cell range
                dims = [1 100]; %input dialog dimensions
                defaults = {sheets{1}, 'All cells'}; %default is sheet 1 and all cells in that sheet
                response = inputdlg(prompt, dlgtitle, dims, defaults, opts);
                if isempty(response) %if user cancels input dialog
                    out_data = false; return %report failure and terminate
                else %if user does not cancel input dialog
                    sheet_name = response{1}; %get sheet name and
                    if strcmp(response{2}, 'All cells'); cell_range = ''; else; cell_range = response{2}; end %cell range
                end
            else %if there is only one sheet in the loaded Excel file
               opts.Interpreter = 'tex';
                opts.Default = 'Do not start Excel';
                response = questdlg(['\fontsize{10}Program will request cell range to read data from. If no range is specified, all data in'...
                    ' the Spreadsheet will be read. Select action.'], ... %inform that program will ask for cell range
                'Reading Excel file','Start Excel', 'Do not start Excel', 'Cancel', opts);
                switch response
                        case 'Start Excel'; open_excel = true; %user pressed Start Excel button
                        case 'Do not start Excel'; open_excel = false; %user pressed Do not start Excel
                        case 'Cancel'; out_data = false; return %user pressed Cancel button
                        case ''; out_data = false;  return %user closed question dialog
                end
                if open_excel; try winopen(file_path); pause(1); catch; open(file_path); pause(1); end; end %try open Excel file in Excel, otherwise in MATLAB
                opts.Interpreter = 'tex';
                prompt = {'\color{blue} \fontsize{10}Enter cell range ( example A1:A40 ). Leave edit field unchanged to read entire sheet.'};
                dlgtitle = 'Specify data location in Excel file'; %ask for cell range
                dims = [1 100];
                defaults = {'All cells'};
                response = inputdlg(prompt, dlgtitle, dims, defaults, opts);
                if isempty(response) %if user closes input dialog
                    out_data = false; return  %terminate execution
                else %otherwise
                    sheet_name = sheets{1};
                    if strcmp(response{1}, 'All cells') %use default cell range if user did not edit field
                        cell_range = '';
                    else
                        try
                            cell_range = response{1}; %otherwise try to extract user input
                        catch
                              opts.Interpreter = 'tex'; opts.WindowStyle = 'modal'; %if that failes, report error
                              errordlg('\fontsize{10}Invalid cell range.', 'Reading Excel file', opts);
                              out_data = false;
                              return
                        end
                    end
                end
            end
            %______________________________________________________________
                out_data = readtable(file_path,'Sheet', sheet_name, 'Range', cell_range, 'ReadVariableNames',false, 'ReadRowNames', false); %read table from file
                out_data = table2array(out_data); %convert table to array
        end
        %------------------------------------------------------------------
        
        %TEXT: txt, DAT data, Comma Separated Values csv
        %------------------------------------------------------------------
        if isfile(file_path) && ( strcmp(file_extension, '.txt') || strcmp(file_extension, '.dat') || strcmp(file_extension, '.csv'))
            format_match = true;%update format match is file extension matches the above
              try
                    out_data = readtable(file_path, 'ReadVariableNames',false, 'ReadRowNames', false); %attempt to read table from file
                    out_data = table2array(out_data); %convert table to array
              catch  %if that fails
                  try
                   fileID = fopen(file_path,'r'); %try to scan the file and read it as 1D data
                   formatSpec = '%f';
                   out_data = fscanf(fileID,formatSpec);
                   fclose(fileID);
                  catch %if that also fails
                       opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
                       errordlg(['\fontsize{10}Failed to read file: ' file_path], 'Reading text file', opts);  %report error
                       out_data = false;
                       return
                  end
              end
        end
        %------------------------------------------------------------------
        
        %MATLAB file: .mat
        if isfile(file_path) && strcmp(file_extension, '.mat')
            format_match = true; %update format match
            try
                out_data = load(file_path, file_name); %assume variable name in mat file is same as file name
                out_data = out_data.file_name;
            catch %if that fails, ask user for variable name to load from mat file
                opts.Interpreter = 'tex';
                prompt = {'\color{blue} \fontsize{10}Enter variable name'};
                dlgtitle = ['Specify variable to load from MAT file: ' file_name];
                dims = [1 100];
                defaults = {''};
                response = inputdlg(prompt, dlgtitle, dims, defaults, opts);
                if isempty(response) %terminate if input dialog is closed
                    out_data = false;
                    return
                else
                    out_data = load(file_path, response{1}); %otherwise load specified variable
                    out_data = out_data.(response{1});
                end            
            end
        end
        
        %Check compliance with expected dimensions
        if ~exist('out_data', 'var'); out_data = false; end %default output
        n_dims = ndims(out_data);
        if n_dims == 2 && min(size(out_data)) == 1; n_dims = 1; end
        if restricted_size && format_match && n_dims ~= expected_dims %if dimensions specified and data does not match
            error_msg = ['\fontsize{10}Expected ' num2str(expected_dims) '-dimensional data. ' ... %raise error
                        'Selected data is ' num2str(n_dims) '-dimensional'];
            opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
            errordlg(error_msg, 'File read error', opts);  
            out_data = false; %alert main program
            return
        end
        
        %Check if data format was recognised
        if ~format_match %if not
            opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
            errordlg('\fontsize{10}File format not recognised.', 'File read error', opts);  %raise error
            out_data = false; %alert main program
            return
        end
        
        %Convert to double
%         if ~isequal(out_data, false); out_data = double(out_data); end
        
end