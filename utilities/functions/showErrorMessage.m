% Displays an error message on the app or prints a warning to the console 
% of the error type and specific error
%
% Copyright (C) 2021 Josiah W. Smith
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
function showErrorMessage(obj,errStr,titleStr)
if nargin == 2
    titleStr = "Warning";
elseif nargin == 1 || nargin == 0
    showErrorMessage(obj,"showErrorMessage must have at least two arguments!","Invalid Function Call");
end
if obj.isApp
    uiconfirm(obj.app.UIFigure,errStr,titleStr,...
        "Options",{'OK'},'Icon','warning');
else
    warning(titleStr + ": " + errStr);
end
end