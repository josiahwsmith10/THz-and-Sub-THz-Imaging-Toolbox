% Checks if the app is using MIMO array or EPC array
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
function tf = isAppEPC(app)
try
    if isApp(app)
        if app.MIMOSwitch.Value == "Use MIMO Array"
            tf = false;
        elseif app.MIMOSwitch.Value == "Use EPC Virtual Elements"
            tf = true;
        end
    else
        tf = false;
    end
catch
    tf = false;
end