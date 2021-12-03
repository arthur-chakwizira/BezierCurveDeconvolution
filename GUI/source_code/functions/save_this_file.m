function save_this_file(save_data)
%        This function accepts the structure save_data which
%        contains the data to save, the folder to save to, the name of the
%        file to save and the format in which to save. This function saves
%        the input data and returns no output.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    data_to_save = save_data.data_to_save;
    save('data_to_save.mat', 'data_to_save')
    this_folder = save_data.this_folder;
    this_name = save_data.this_name;
    this_format = save_data.this_format;
    %Check if input data is 1-D and make sure it will be saved as a column
    if ismatrix(data_to_save) && size(data_to_save, 1) == 1; data_to_save = data_to_save'; end
    %----------------------------------------------------------------------
    if isequal(this_format, 0) && isequal(this_name,0) %this will happen when the result is from uiputfile
        [path_to_folder, file_name, file_format] = fileparts(this_folder);
        if isempty(file_format) %may happen if user chooses to save in format "All files"
            file_format = '.nii';
            opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
            msgbox('\fontsize{10}No format provided for data saving. Data has been saved to format NIFTI (*.nii)', 'File save info','warn' , opts); 
        end
    else %when function is called from script; all arguments provided
        path_to_folder = this_folder;
        file_name = this_name;
        switch this_format
            case 'NIFTI (*.nii)'; file_format = '.nii';
            case 'DICOM (*.dcm)'; file_format = '.dcm'; %maybe later
            case 'MATLAB (*.mat)'; file_format = '.mat';
            case 'Text (*.txt)'; file_format = '.txt';
            case 'DAT (*.dat)'; file_format = '.dat';
            case 'Comma Separated Values (*.csv)'; file_format = '.csv';
            case 'Excel (*.xlsx)'; file_format = '.xlsx';
        end
    end
    
    %----------------------------------------------------------------------
    switch file_format
        case '.nii'; niftiwrite(data_to_save, strcat(path_to_folder,'\', file_name, file_format));
        case '.dcm'
            header_info = save_data.header_info;
            size_data_to_save = size(data_to_save);
            switch ndims(data_to_save)
                case 2 %if data has is mxn or 1xn
                    if ~isempty(size_data_to_save(size_data_to_save == 1))  %if it is 1xn
                        opts.Interpreter = 'tex'; opts.WindowStyle = 'modal';
                        errordlg('\fontsize{10}Can not save 1-D data to format DICOM (*.dcm)', 'File save error', opts); %raise this error 
                    else %if it is mxn
                        save_destination = strcat(path_to_folder, '\', file_name, file_format);
                        header_info.Filename = save_destination; %update header file from loaded data 
                        try
                            dicomwrite(data_to_save, save_destination, header_info); %try to save dicom file with that header file
                        catch %if header verification fails
                            try
                                dicomwrite(data_to_save, save_destination, header_info, 'CreateMode', 'copy'); % disregard the verification
                            catch %if even that fails
                                dicomwrite(data_to_save, save_destination) %let matlab generate generic header file
                            end
                        end
                    end
                case 3
                    n_slices = size_data_to_save(3);
                    for sl = 1:n_slices
                        save_destination = strcat(path_to_folder, '\', file_name,'_', num2str(sl));
                        header_info.Filename = save_destination; %update header file from loaded data 
                        try
                            dicomwrite(squeeze(data_to_save(:,:,sl)), save_destination, header_info)
                        catch
                            try
                                dicomwrite(squeeze(data_to_save(:,:,sl)), save_destination, header_info, 'CreateMode', 'copy')
                            catch
                                dicomwrite(squeeze(data_to_save(:,:,sl)), save_destination)
                            end
                        end
                    end
                case 4
                    save_destination = strcat(path_to_folder, '\', file_name, file_format);
                    header_info.Filename = save_destination; %update header file from loaded data 
                    try
                        dicomwrite(data_to_save, save_destination, header_info)
                    catch
                        try
                            dicomwrite(data_to_save, save_destination, header_info, 'CreateMode', 'copy')
                        catch
                            dicomwrite(data_to_save, save_destination)
                        end
                    end
            end
            
        case '.mat'; save(strcat(path_to_folder, '\', file_name, file_format), data_to_save);
        case '.txt'
%             fileID = fopen(strcat(path_to_folder, '\', file_name, file_format), 'w');
%             fprintf(fileID,'%f\n', data_to_save);
%             fclose(fileID);
               if ~ismatrix(data_to_save)
                  opts.Interpreter = 'tex';
                  opts.Default = 'Save to NIFTI';
                  response = questdlg(['\fontsize{10}Attempting to save ' num2str(ndims(data_to_save)) ' - dimensional data to format: Text (*.txt)'...
                    ' This may result in generation of a large number of files. Select action.'], ...
                    'Saving to text file','Continue', 'Save to NIFTI', 'Save to DICOM', opts);
                  switch response
                      case 'Continue'
                          switch ndims(data_to_save)
                              case 3
                                  for sl = 1:size(data_to_save, 3)
                                      slice_to_save = array2table(squeeze(data_to_save(:,:,sl)));
                                      writetable(slice_to_save, strcat(path_to_folder, '\', file_name, '_Slice_',num2str(sl), file_format), 'Delimiter', 'tab', 'WriteVariableNames', false,'WriteRowNames', false )
                                  end
                              case 4
                                  for sl = 1:size(data_to_save, 3)
                                      for tp = 1:size(data_to_save, 4)
                                          array_to_save = array2table(squeeze(data_to_save(:,:,sl, tp)));
                                          writetable(array_to_save, strcat(path_to_folder, '\', file_name, '_Slice_',num2str(sl), '_Image_', num2str(tp), file_format), 'Delimiter', 'tab', 'WriteVariableNames', false,'WriteRowNames', false )
                                      end
                                  end
                          end
                      case 'Save to NIFTI' %HERE HERE
                          save_data.this_format = 'NIFTI (*.nii)';
                            save_this_file(save_data)
                      case 'Save to DICOM'
                          save_data.this_format = 'DICOM (*.dcm)';
                            save_this_file(data_to_save, path_to_folder, file_name, 'DICOM (*.dcm)')
                      case ''
                          return
                  end  
               else 
                data_to_save = array2table(data_to_save);
                writetable(data_to_save, strcat(path_to_folder, '\', file_name, file_format), 'Delimiter', 'tab', 'WriteVariableNames', false,'WriteRowNames', false )
               end
        case {'.dat' ,'.csv'}
                  if ~ismatrix(data_to_save)
                  opts.Interpreter = 'tex';
                  opts.Default = 'Save to NIFTI';
                  response = questdlg(['\fontsize{10}Attempting to save ' num2str(ndims(data_to_save)) ' - dimensional data to format: ' this_format ...
                    ' This may result in generation of a large number of files. Select action.'], ...
                    ['Saving to' this_format 'file'],'Continue', 'Save to NIFTI', 'Save to DICOM','Cancel', opts);
                  switch response
                      case 'Continue'
                          switch ndims(data_to_save)
                              case 3
                                  for sl = 1:size(data_to_save, 3)
                                      slice_to_save = array2table(squeeze(data_to_save(:,:,sl)));
                                      writetable(slice_to_save, strcat(path_to_folder, '\', file_name, '_Slice_',num2str(sl), file_format), 'Delimiter', 'tab', 'WriteVariableNames', false,'WriteRowNames', false )
                                  end
                              case 4
                                  for sl = 1:size(data_to_save, 3)
                                      for tp = 1:size(data_to_save, 4)
                                          array_to_save = array2table(squeeze(data_to_save(:,:,sl, tp)));
                                          writetable(array_to_save, strcat(path_to_folder, '\', file_name, '_Slice_',num2str(sl), '_Image_', num2str(tp), file_format), 'Delimiter', 'tab', 'WriteVariableNames', false,'WriteRowNames', false )
                                      end
                                  end
                          end
                      case 'Save to NIFTI' %HERE HERE
                            save_this_file(data_to_save, path_to_folder, file_name, 'NIFTI (*.nii)')
                      case 'Save to DICOM'
                            save_this_file(data_to_save, path_to_folder, file_name, 'DICOM (*.dcm)')
                      case 'Cancel'
                          return
                  end  
                  else    
                    data_to_save = array2table(data_to_save);
                    writetable(data_to_save, strcat(path_to_folder, '\', file_name, file_format), 'WriteVariableNames', false)
                  end
            
        case {'.xlsx' , '.xls'}
            if ~ismatrix(data_to_save)
                          switch ndims(data_to_save)
                              case 3
                                  for sl = 1:size(data_to_save, 3)
                                      slice_to_save = squeeze(data_to_save(:,:,sl));
                                      write_this_matrix(slice_to_save, strcat(path_to_folder, '\', file_name, file_format), 'Sheet', ['Slice_' num2str(sl)])
                                  end
                              case 4
                                  for sl = 1:size(data_to_save, 3)
                                      for tp = 1:size(data_to_save, 4)
                                          array_to_save = squeeze(squeeze(data_to_save(:,:,sl, tp)));
                                          write_this_matrix(array_to_save, strcat(path_to_folder, '\', file_name,'_Slice_',num2str(sl), file_format), 'Sheet', ['Image_' num2str(tp)])
                                      end
                                  end
                          end                
            else
            write_this_matrix(data_to_save, strcat(path_to_folder, '\', file_name, file_format))
            end
    end
 
end