function header = get_dicom_header(path_to_data)
%        This function accepts a path to a dicom file/series and returns the header
%        information in the file/series.
%        Author: 
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if isfile(path_to_data) %means user chose dicom 4d volume
        else
        if isfolder(path_to_data) %if it is a folder
            folder_contents = dir(path_to_data); %get all files in the folder
            file_names = {folder_contents.name}; %get file names of these files
            file_to_use = fullfile(path_to_data, file_names{round(length(file_names)/2)}); %choose the middle file in the series
            path_to_data = file_to_use; %and extract its part
        end
  end
    try
       header = dicominfo(path_to_data); %attempt to read its header
    catch
       header = false; %otherwise return false
    end
end