function [unique_fname, change_detect] = auto_rename (full_file_str, appendage_str)

%unique_fname = auto_rename (full_file_str, appendage_str)
%   auto_rename will automatically rename any file if the filename already exists in the directory chosen. The renamed file
%    will contain a numerical tag as described below.  The user may use the output to name and save the file.  
%    Can be used to name mat files, figures, m files, any file.  
%
%  INPUTS:    
%   full_file_str is a string containing the path & filename & extension  Example: 'C:\Users\adanz\Documents\MATLAB\filename.mat'
%
%   appendage_str is a string that describes how duplicate files should be appended.  Use the number '0' to indicate 
%    where the numerical tag should go.  The numerical tag will start from 2 and will search for the next available unique filename. 
%       Examples: ' (0)'        will tag duplicate file, 'filename.mat' with 'filename (2).mat' and then filename (3).mat etc.....
%                 '_0'          will tag duplicate file, 'filename.mat' with 'filename_2.mat'
%                 ' copy (0)'   will tag duplicate file, 'filename.mat' with 'filename copy (2).mat'
%                 '0'           will tag duplicate file, 'filename.mat' with 'filename2.mat
%       
%  OUTPUTS:   
%   unique_fname is a full path/filename/ext in str format.
%   change_detect is logical (1/0) declaring if the original full_file_str was changed due to duplication (1=changed)
%
%  EXAMPLES:
%   unique_fname = auto_rename('C:\Users\adanz\Documents\MATLAB\solar_data.mat', '_(0)')
%                = 'C:\Users\adanz\Documents\MATLAB\solar_data_(2).mat'
%                = 'C:\Users\adanz\Documents\MATLAB\solar_data_(3).mat'  etc...
%
%   unique_fname = auto_rename('C:\Users\adanz\Documents\MATLAB\lunar_dist.fig', ' vs.0')
%                = 'C:\Users\adanz\Documents\MATLAB\lunar_dist vs2.fig'
%            ... = 'C:\Users\adanz\Documents\MATLAB\lunar_dist vs9.fig'
%
%   To parse the output into path/file name/extension, use fileparts(unique_fname).
%   To save file, use save(unique_fname).
%
%  Adam Danz 140412
%%
[PATHSTR,NAME,EXT] = fileparts(full_file_str);
a = 2;                                              %'a' is the file appendate number
z = 0;                                              %'z' is '0' marker 
change_detect = logical(0);

%make sure dir exists
if logical(exist(PATHSTR, 'dir') ~= 7)
    disp ([PATHSTR, ' does not exist'])
    return
end

%make sure ext is present
if isempty(EXT)
    disp ('File extension must be present in auto-rename inputs.')
    return
end

%make sure NAME is present
if isempty(NAME)
    disp ('File name must be present in auto-rename inputs along with file extension.')
    return
end

%make sure appendage_str is present
if nargin <2 || isempty(appendage_str)
    disp ('Appendage string missing.  Examples:  ''_0'', ''vs.0'', '' (0)'', etc.')
    return
end

%make sure appendage contains '0'
if isempty(strfind(appendage_str, '0'))
    disp ('Appendage string missing ''0'' which indicates where the numerical appendage should go.')
    disp ('Examples:  ''_0'', ''vs.0'', '' (0)'', etc.')
    return
end

%make sure appendage doesn't contain period which would screw up the extension (replaces '.' with []
appendage_str(strfind(appendage_str, '.')) = [];

%Check for existence of file and tag on appendage if necessary
tempname = NAME;
while exist(full_file_str,'file') == 2
    appendage_str = strrep(appendage_str, num2str(z), num2str(a));
    tempNAME = [NAME, appendage_str];
    full_file_str = fullfile(PATHSTR, [tempNAME, EXT]);
    z = a;
    a = a+1;
    change_detect = logical(1);
end

%output
unique_fname = full_file_str;
    
%% Notes

% % use this code to create a bunch of mat files
% f = 'C:\Users\adanz\Documents\MATLAB\savehere\dragons den\filename.mat'
% f2 = auto_rename (f, ' (0)')
% save (f2)


