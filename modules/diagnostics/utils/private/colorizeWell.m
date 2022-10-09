function c = colorizeWell(type, index, D)
%Undocumented Utility Function

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    switch lower(type)
        case 'prod'
            if numel(D.prod)<6
                cmap = lines(numel(D.prod));
            else
                cmap = jet(numel(D.prod));
            end
            c = cmap(index,:);
        case 'inj'
            cmap = gray(numel(D.inj));
            c = cmap(index,:);
        case 'global'
            if any(D.prod == index)
                c = colorizeWell('prod', find(D.prod == index), D);
            else
                c = colorizeWell('inj', find(D.inj == index), D);
            end
    end
end
