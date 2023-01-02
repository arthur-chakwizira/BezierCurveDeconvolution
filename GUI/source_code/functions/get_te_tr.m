function [te, tr] = get_te_tr(handles)
%       This function accepts the handles structure and returns Echo time
%       and Repetition time
%        Author:
%              Arthur Chakwizira
%              arthurchakwizira@gmail.com
%             Medical Radiation Physics, Lund University, Sweden
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_folder = handles.file_folder;
te_tr_fn = fullfile(file_folder, "te_tr.txt");
if isfile(te_tr_fn)
    fileID = fopen(te_tr_fn,'r'); %try to scan the file and read it as 1D data
    formatSpec = '%f';
    te_tr = fscanf(fileID,formatSpec);
    fclose(fileID);
    if ~isempty(te_tr)
    te = te_tr(1);
    tr = te_tr(2);
    else
            te = [];
    tr = [];
    end
else
    te = [];
    tr = [];
end

wrap_text = 'Waiting for user input...';
set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'b')
opts.Interpreter = 'tex';
prompt = {'\color{blue} \fontsize{10} TE [ms] :', '\color{blue} \fontsize{10} TR [s] :'};
dlgtitle = 'DECONVOLVER: Confirm TE/TR';
dims = [1 80];

current_te = num2str(te*1E3);
current_tr = num2str(tr);

defaults = {current_te, current_tr};

response = inputdlg(prompt, dlgtitle, dims, defaults, opts);

if isempty(response) %if use clicks cancel, report and terminate
    wrap_text = 'Analysis terminated by user. Re-run to select data.';
    set(handles.edit7, 'String', wrap_text, 'ForegroundColor', 'r')
    te = false; tr = false;
    return
end

te = str2double(response{1})*1E-3;
tr = str2double(response{2});

fileID = fopen(te_tr_fn,'w'); %
formatSpec = '%4.4f\n';
fprintf(fileID,formatSpec, te);
fprintf(fileID,formatSpec, tr);
fclose(fileID);

end
