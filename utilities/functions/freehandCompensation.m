% Computes the freehand MIMO compensation derived by Josiah
%
% Copyright (C) 2021 Josiah W. Smith
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.

function obj = freehandCompensation(obj)
% Y-dimension only multistatic-to-monostatic approximation
if numel(obj.scanner.sarSize) == 3
    k = reshape(obj.k_vec,1,1,[]);
elseif numel(obj.scanner.sarSize) == 4
    k = reshape(obj.k_vec,1,1,1,[]);
else
    warning("Could not perform multistatic-to-monostatic approximation!!")
    obj.sarData = single(zeros(size(reshape(obj.sarData,[],size(obj.sarData,3),size(obj.sarData,4)))));
    return;
end

obj.sarData = squeeze(reshape(obj.sarData,obj.ant.vx.numVx,obj.scanner.numY,obj.scanner.numX,obj.wav.Nk));

dly = obj.ant.vx.dxy(:,2);

dlz = -reshape(obj.scanner.vx.xyz_m(:,3),obj.ant.vx.numVx,obj.scanner.numY,obj.scanner.numX) + obj.ant.z0_m;

phi_l = 2*dlz + (dly.^2)/(4*obj.zRef_m);

compensationTerm = exp(-1j*k.*phi_l);

obj.sarData = obj.sarData .* squeeze(compensationTerm);
obj.sarData(isnan(obj.sarData)) = 0;

obj.sarData = squeeze(reshape(obj.sarData,obj.ant.vx.numVx*obj.scanner.numY,obj.scanner.numX,obj.wav.Nk));

end