function target_folder = generate_target_folder(handles)
c = string(datetime());
c = strrep(c, ":", "_");
c = strrep(c, " ", " at ");
target_folder = fullfile(handles.file_folder, handles.file_name, c); %save in folder with same name asa input file
if ~exist(target_folder, 'dir'); mkdir(target_folder); end %if target folder is non-existent, make it
end