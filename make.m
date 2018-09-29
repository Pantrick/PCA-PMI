function make(varargin)
%if you want to remove the path from matlab search path, go to 'Set Path'
% Example:
% `make` 
% `make -p` Path will be added to Matlab search path permenently
    permenent = 0;
    if any(cellfun(@(s)strcmp(s,'-p'),varargin))
        disp('Path will be added to Matlab search path permenently');
        disp('If you want to remove the path from matlab search path, please go to ''Set Path''');
        permenent = 1;
    end
    addpath(genpath('lib'));
    if permenent
        savepath;
    end
end