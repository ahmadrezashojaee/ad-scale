function x = recoverVars(eq, n, sol)
% Recover previously eliminated variables x at position n using solutions sol

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
   solInx = [1:(n-1), (n+1):(numel(sol)+1)];
   x = - eq.jac{n}\(eq.val);
   for k  = 1:numel(solInx)
       if(~isempty(x) && ~isempty(sol{k}))
           x = x - eq.jac{n}\(eq.jac{solInx(k)}*sol{k});
       end
   end
end
