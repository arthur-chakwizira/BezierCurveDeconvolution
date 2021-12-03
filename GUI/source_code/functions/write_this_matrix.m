function writematrix(A, filename, varargin)
%WRITEMATRIX Write a matrix to a file.
%   WRITEMATRIX(A) writes the homogenous array A to a comma-delimited text
%   file. The file name is the workspace name of the homogenous array A,
%   appended with '.txt'. If WRITEMATRIX cannot construct the file name
%   from the homogenous array input, it writes to the file 'matrix.txt'.
%   WRITEMATRIX overwrites any existing file.
%
%   WRITEMATRIX(A,FILENAME) writes the homogenous array A to the file
%   FILENAME as column-oriented data. WRITEMATRIX determines the file
%   format from its extension. The extension must be one of those listed
%   below.
%
%   WRITEMATRIX(A,FILENAME,'FileType',FILETYPE) specifies the file type,
%   where FILETYPE is one of 'text' or 'spreadsheet'.
%
%   WRITEMATRIX writes data to different file types as follows:
%
%   .txt, .dat, .csv:  Delimited text file (comma-delimited by default).
%
%          WRITEMATRIX creates a column-oriented text file, i.e., each
%          column of each variable in A is written out as a column in the
%          file.
%
%          Use the following optional parameter name/value pairs to control
%          how data are written to a delimited text file:
%
%          'Delimiter'      The delimiter used in the file. Can be any of ' ',
%                           '\t', ',', ';', '|' or their corresponding names 'space',
%                           'tab', 'comma', 'semi', or 'bar'. Default is ','.
%
%          'QuoteStrings'   A logical value that specifies whether to write
%                           text out enclosed in double quotes ("..."). If
%                           'QuoteStrings' is true, any double quote characters that
%                           appear as part of a text variable are replaced by two
%                           double quote characters.
%
%          'DateLocale'     The locale that writematrix uses to create month and
%                           day names when writing datetimes to the file. LOCALE must
%                           be a character vector or scalar string in the form xx_YY.
%                           See the documentation for DATETIME for more information.
%
%          'Encoding'       The encoding to use when creating the file.
%                           Default is 'system' which means use the system's default
%                           file encoding.
%
%   .xls, .xlsx, .xlsb, .xlsm, .xltx, .xltm:  Spreadsheet file.
%
%          WRITEMATRIX creates a column-oriented spreadsheet file, i.e., each column
%          of each variable in A is written out as a column in the file.
%
%          Use the following optional parameter name/value pairs to control how data
%          are written to a spreadsheet file:
%
%          'DateLocale'     The locale that writematrix uses to create month and day
%                           names when writing datetimes to the file. LOCALE must be
%                           a character vector or scalar string in the form xx_YY.
%                           Note: The 'DateLocale' parameter value is ignored
%                           whenever dates can be written as Excel-formatted dates.
%
%          'Sheet'          The sheet to write, specified the worksheet name, or a
%                           positive integer indicating the worksheet index.
%
%          'Range'          A character vector or scalar string that specifies a
%                           rectangular portion of the worksheet to write, using the
%                           Excel A1 reference style.
%
%          'UseExcel'      A logical value that specifies whether or not to create the
%                          spreadsheet file using Microsoft(R) Excel(R) for Windows(R). Set
%                          'UseExcel' to one of these values:
%                               true  -  Opens an instance of Microsoft
%                                        Excel to write (or read) the file.
%                                        This is the default setting for
%                                        Windows systems with Excel
%                                        installed.
%                               false -  Does not open an instance of
%                                        Microsoft Excel to write (or read)
%                                        the file. Using this setting may
%                                        cause the data to be written
%                                        differently for files with live
%                                        updates (e.g. formula evaluation
%                                        or plugins).
%
%   In some cases, WRITEMATRIX creates a file that does not represent A
%   exactly, as described below. If you use READMATRIX(FILENAME) to read that
%   file back in and create a new matrix, the result may not have exactly
%   the same format or contents as the original matrix.
%
%   *  WRITEMATRIX writes out numeric data using long g format, and
%      categorical or character data as unquoted text.
%   *  WRITEMATRIX writes out arrays that have more than two dimensions as two
%      dimensional arrays, with trailing dimensions collapsed.
%
%   See also READMATRIX, READTABLE, READCELL, WRITETABLE, WRITECELL.

%   Copyright 2018 The MathWorks, Inc.


matrixTypes = ["duration","datetime","categorical","string","char","logical"];
if ~isnumeric(A) && ~any(strcmp(class(A),matrixTypes))
    if iscell(A)
        error(message('MATLAB:table:write:UnsupportedInputType','cell','writecell'));
    elseif istable(A)
        error(message('MATLAB:table:write:UnsupportedInputType','table','writetable'));
    elseif istimetable(A)
        error(message('MATLAB:table:write:UnsupportedInputType','timetable','writetimetable'));
    else
        error(message('MATLAB:table:write:UnsupportedTypeIn',class(A)));
    end
end

if nargin < 2
    matrixname = inputname(1);
    if isempty(matrixname)
        matrixname = 'matrix';
    end
    filename = [matrixname '.txt'];
else
    for i = 1:2:numel(varargin)
        n = strlength(varargin{i});
        if n > 5 && strncmpi(varargin{i},'WriteVariableNames',n)
            error(message('MATLAB:table:write:WriteVariableNamesNotSupported','WRITEMATRIX'));
        end
        if n > 5 && strncmpi(varargin{i},'WriteRowNames',n)
            error(message('MATLAB:table:write:WriteRowNamesNotSupported','WRITEMATRIX'));
        end
    end
end

try
    T = table(A);
    writetable(T,filename,varargin{:},'WriteVariableNames',false);
catch ME
    throw(ME)
end



end
