function [fn,pn] = myuiputfile(ext,title)

%****************************************
% Part of NIRFAST
% Last Updated 8/10/09 - M Jermyn
%****************************************
%
% [fn,pn] = myuiputfile(ext,title)
%
% Allows the user to browse for saving, remembers last path.
%
% ext is the file extension
% title is the title of the window
% fn is the resulting filename
% pn is the resulting path


% retrieve last used path
if ispref('nirfast_gui','lastpath')
    lastpath = getpref('nirfast_gui','lastpath');
else
    lastpath = '/';
end

% get file
[fn,pn] = uiputfile(ext,title,lastpath);

% save path
if pn
    setpref('nirfast_gui','lastpath',pn);
end