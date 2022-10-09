function bi = blockInverter(opt)
%Undocumented Utility Function

%{
Copyright 2020 University of Bergen and SINTEF Digital, Mathematics & Cybernetics.

This file is part of the MPSA-W module for the MATLAB Reservoir Simulation Toolbox (MRST).

The MPSA-W module is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

The MPSA-W module is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with the MPSA-W module.  If not, see <http://www.gnu.org/licenses/>.
%}


   if ~ischar(opt.invertBlocks)
      dispif(opt.verbose, ...
            ['Unsupported option value of type ''%s'' in ', ...
             'option ''invertBlocks''. Reset to default ' , ...
             '(''matlab'')\n'], class(opt.invertBlocks));
      opt.invertBlocks = 'matlab';
   end

   switch lower(opt.invertBlocks)
      case {'matlab', 'm', 'builtin'}
         bi = @invertDiagonalBlocks;
      case {'mex', 'c', 'accelerated'}
         bi = @invertDiagonalBlocksMex;
      otherwise
         dispif(opt.verbose, ...
               ['Unsupported value ''%s'' in option ', ...
                '''invertBlocks''.\nMust be one of ''matlab'' or ', ...
                '''mex''.\nReset to default (''matlab'').'], ...
                opt.invertBlocks);

         bi = @invertDiagonalBlocks;
   end
end

% Matlab code calling mex function invertSmallMatrices.c
function iA = invertDiagonalBlocksMex(A, sz)
   sz     = int32(sz);
   blocks = matrixBlocksFromSparse(A, sz);
   iA     = blockDiagMatrix(invv(blocks, sz), sz);
end

%--------------------------------------------------------------------------

% Pure Matlab code using inv
function iA = invertDiagonalBlocks(A, sz)
   V = zeros([sum(sz .^ 2), 1]);
   [p1, p2] = deal(0);

   for b = 1 : numel(sz)
      n  = sz(b);
      n2 = n * n;
      i  = p1 + (1 : n);

      V(p2 + (1 : n2)) = inv(full(A(i, i)));

      p1 = p1 + n;
      p2 = p2 + n2;
   end

   iA = blockDiagMatrix(V, sz);
end

%--------------------------------------------------------------------------

function A = blockDiagMatrix(V, sz)
   [I, J] = blockDiagIndex(sz);
   A      = sparse(I, J, V);
end
