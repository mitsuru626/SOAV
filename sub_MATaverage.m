function DataMean = sub_MATaverage(target_regexp, varname, varargin)
% sub_MATaverage - Compute average of variable in Mat-files
%
% DataMean = sub_MATaverage(target_regexp,varname,option)
%
% Arguments:
%   target_regexp: regular expression of target files
%   varname: name of target variable
% Option:
%   SaveVarName: name of variable assigned in workspace
%
% Mitsuru Toyoda (Tokyo Metropolitan University)
% 

DataMean = 0;
n_data = 0;
data_list = dir(target_regexp);
if isempty(data_list)
    error('No file found !');
end
fprintf(['File: "', target_regexp, '" / Varname: "', varname, '"\n']);
for i = 1:length(data_list)
    % output of load command is structure
    if isfield(load(data_list(i).name, '-mat'), varname)
        data_i = getfield(load(data_list(i).name, '-mat', varname), varname);
        fprintf([num2str(i), '/', num2str(length(data_list)), ': ', data_list(i).name, '\n']);
        DataMean = DataMean + data_i;
        n_data = n_data + 1;
    end
end
DataMean = DataMean / n_data;

% assign variable in workspace
if ~isempty(varargin) && n_data ~= 0
    assignin('base', varargin{1}, DataMean);
end

end