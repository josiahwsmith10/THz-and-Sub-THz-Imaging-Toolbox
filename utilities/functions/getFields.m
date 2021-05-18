% Sets properties of savedobj to obj and returns obj
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
function obj = getFields(obj,savedobj,skipfields)
% Inputs
%   obj         -   Object to set parameters to
%   savedobj    -   Object to get parameters from
%   skipfields  -   1D string array of parameter field names to skip

if nargin == 2
    skipfields = [];
end
savedfields = string(fieldnames(savedobj));
for indField = 1:length(savedfields)
    currfield = savedfields(indField);
    if max(currfield == string(skipfields(:)))
        continue
    else
        obj.(currfield) = savedobj.(currfield);
    end
end
end