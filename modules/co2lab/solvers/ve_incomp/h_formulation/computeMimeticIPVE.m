function S = computeMimeticIPVE(G, rock, varargin)
%Compute mimetic inner product matrices.
%
% SYNOPSIS:
%   S = computeMimeticIPVE(G, rock)
%   S = computeMimeticIPVE(G, rock, 'pn', pv, ...)
%
% PARAMETERS:
%   G       - Grid structure made from function 'topSurfaceGrid'.
%
%   rock    - Rock data structure with valid field 'perm'.  The
%             permeability is assumed to be in measured in units of
%             metres squared (m^2).  Use function 'darcy' to convert from
%             (milli)darcies to m^2, e.g.,
%
%                 perm = convertFrom(perm, milli*darcy)
%
%             if the permeability is provided in units of millidarcies.
%
%             The field rock.perm may have ONE column for a scalar
%             permeability in each cell, TWO/THREE columns for a diagonal
%             permeability in each cell (in 2/3 D) and THREE/SIX columns
%             for a symmetric full tensor permeability.  In the latter
%             case, each cell gets the permeability tensor
%
%                 K_i = [ k1  k2 ]      in two space dimensions
%                       [ k2  k3 ]
%
%                 K_i = [ k1  k2  k3 ]  in three space dimensions
%                       [ k2  k4  k5 ]
%                       [ k3  k5  k6 ]
%
%   'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%             supported options are:
%               - Type -- The kind of system to assemble.  The choice made
%                         for this option influences which pressure solvers
%                         can be employed later.
%                         String.  Default value = 'hybrid'.
%                 Supported values are:
%                   - 'hybrid'      : Hybrid system with inverse of B
%                                     (for Schur complement reduction)
%                   - 'mixed'       : Hybrid system with B
%                                     (for reduction to mixed form)
%                   - 'tpfa'        : Hybrid system with B
%                                     (for reduction to tpfa form)
%                   - 'comp_hybrid' : Both 'hybrid' and 'mixed'
%
%               - InnerProduct -- The inner product with which to define
%                                 the mass matrix.
%                                 String.  Default value = 'ip_simple'.
%                 Supported values are:
%                   - 'ip_simple'   : inner product having the 2*tr(K)-term.
%                   - 'ip_tpf'      : inner product giving method equivalent
%                                     to two-point flux approximation (TPFA).
%                   - 'ip_quasitpf' : inner product ''close to'' TPFA
%                                     (equivalent for Cartesian grids with
%                                      diagonal permeability tensor).
%                   - 'ip_rt'       : Raviart-Thomas for Cartesian grids,
%                                     else not valid.
%
%               - Verbose -- Whether or not to emit progress reports during
%                            the assembly process.
%                            Logical.  Default value = FALSE.
%
%               - facetrans -- If facetrans is specified, the innerproduct
%                              is modified to account for a face
%                              transmissibilities specified as [faces,
%                              face_transmissibilities]
%
%               - intVert  -- Whether or not to compute mobilites
%                             by vertically integrating 3D model (instead
%                             of using averaged values). If false, use
%                             average peremability value of column
%                             (rock2D). Default value: intVert = false.
%
% RETURNS:
%   S - Pressure linear system structure having the following fields:
%         - BI / B : inverse of B / B in hybrid/mixed system - type   :
%         system type (hybrid or mixed)
%         - ip     : inner product name
%
% COMMENTS:
%  In the hybrid discretization, the matrices B, C and D appear as
%       [  B   C  D  ]
%       [  C'  O  O  ]
%       [  D'  O  O  ]
%
% SEE ALSO:
%   `solveIncompFlowVE`, `darcy`, `permTensor`.

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

mrstNargInCheck(2, [], nargin);

opt = struct('InnerProduct', 'ip_simple', ...
             'verbose',      false,       ...
             'Type',         'hybrid',    ...
             'facetrans',    zeros(0, 2), ...
             'intVert',      false,       ...
             'inner_param',  6);
opt = merge_options(opt, varargin{:});

dim = size(G.nodes.coords, 2);
%--------------------------------------------------------------------------
%- Model and discretisation parameters ------------------------------------
%
perm = permTensor(rock, dim);

assert (size(perm,1) == G.cells.num, ...
       ['Permeability must be defined in active cells only.\n', ...
        'Got %d tensors, expected %d (== number of cells).'],   ...
        size(perm,1), G.cells.num);

if opt.intVert % permeability will be included in pseudo mobilites
 perm(:) = 1;
end

[systemType, ip, ipname, computeIP, computeInverseIP] = setup(opt);

%--------------------------------------------------------------------------
%- Compute discrete gradient (C) and flux continuity (D) matrices ---------
%
cf = double(G.cells.faces(:,1));

if opt.verbose,
   fprintf('Computing cell inner products ...\t\t')
   tic
end

cellno = rldecode(1 : G.cells.num, diff(G.cells.facePos), 2) .';
sgn    = 2*double(G.faces.neighbors(cf, 1) == cellno) - 1;

%--------------------------------------------------------------------------
%- Derive geometric primitives per half-contact in given model ------------
%
nc = [G.nodes.coords G.nodes.z];
fc = [G.faces.centroids G.faces.z];
cc = [G.cells.centroids G.cells.z];
e_ij = nc(G.faces.nodes(1:2:end-1),:) - nc(G.faces.nodes(2:2:end),:);
areas        = sqrt(sum(e_ij.^2,2));
areas        = areas(cf);
centVecs_tmp = fc(cf,:) - cc(cellno,:);
surfBasis    = cross(rldecode(G.cells.normals, diff(G.cells.facePos)), ...
                     centVecs_tmp);
cfaceNormals = zeros(size(surfBasis,1),2);
centVecs     = zeros(size(surfBasis,1),2);

for i=1:G.cells.num
   cfind = G.cells.facePos(i):G.cells.facePos(i+1)-1;
   C3D   = centVecs_tmp(cfind,:);
   ee3D  = e_ij(cf(cfind),:);
   U     = orth(surfBasis(cfind,:)')';
   if (dot(cross(U(1,:),U(2,:)), G.cells.normals(i,:)) < 0)
      U = U([2 1],:);
      assert(dot(cross(U(1,:),U(2,:)),G.cells.normals(i,:))>0);
   end
   C2D  = U*C3D';
   ee2D = (U*ee3D')';
   N2D  = [-ee2D(:,2),ee2D(:,1)];
   centVecs(cfind,:) = C2D';
   cfaceNormals(cfind,:) = N2D;
end

outNormals = bsxfun(@times, cfaceNormals, sgn);

clear fc nc cc cf cellno sgn

%--------------------------------------------------------------------------
%- Preallocate SPARSE storage for BI, one value per half-contact ----------
%
dimProd = double(diff(G.cells.facePos));
sizVec  = sum(dimProd .^ 2);

if computeIP,        val_ip    = zeros([sizVec, 1]); end
if computeInverseIP, val_invip = zeros([sizVec, 1]); end

%--------------------------------------------------------------------------
%- Compute half-contact numbering for each cell (-> SPARSE assembly) ------
%

% Current index to row/column in BI (or B), and current index to value.
ix = 0; ix2 = 0;

% Compute half-face transmissibilities from face transmissibilities.
if any(opt.facetrans),
   if ~all(opt.facetrans(:,2)>0),
      error('All face transmissibilities must be positive.');
   end
   hft = zeros(G.faces.num, 1);
   hft(opt.facetrans(:,1)) = opt.facetrans(:,2);
   hft = hft.*sum(G.faces.neighbors ~= 0, 2);
   i=hft~=0; hft(i) = 1./hft(i);
   hft = hft(G.cells.faces(:,1));
else
   hft = zeros(size(G.cells.faces(:,1)));
end

%--------------------------------------------------------------------------
% Main loop -- Compute BI and/or B for each cell in model -----------------
%

for k = 1 : G.cells.num,
   nF  = dimProd(k);
   pR  = ix  + 1 : ix  + nF;   % half-face indices
   pR2 = ix2 + 1 : ix2 + nF^2;

   a = areas(pR);
   v = G.cells.volumes(k);

   K = reshape(perm(k,:), [dim, dim]);
   C = centVecs  (pR,:);
   N = outNormals(pR,:);

   % Compute inner product if caller specified 'mixed' method.
   if computeIP,
      val_ip(pR2)    = ip(a, v, K, C, N, false) + diag(hft(pR));
   end

   % Compute inverse inner product if caller specified 'hybrid' method.
   if computeInverseIP,
      this_ip = ip(a, v, K, C, N, true);
      if any(hft(pR)),
         val_invip(pR2) = inv( inv(this_ip) + diag(hft(pR)) );
      else
         val_invip(pR2) = this_ip;
      end
   end

   ix  = ix  + nF;
   ix2 = ix2 + nF.^2;
end
tocif(opt.verbose)

%--------------------------------------------------------------------------
%- Assemble final block diagonal BI and/or B for complete model -----------
%
[ind1, ind2] = blockDiagIndex(dimProd, dimProd);

clear areas dimProd centVecs outNormals perm
if opt.verbose,
   fprintf('Assembling global inner product matrix ...\t')
   tic
end
if computeIP,
   n    = size(G.cells.faces, 1);
   S.B  = sparse(ind1, ind2, val_ip,    n, n);
end
if computeInverseIP,
   n    = size(G.cells.faces, 1);
   S.BI = sparse(ind1, ind2, val_invip, n, n);
end
tocif(opt.verbose)

if opt.verbose && computeIP && computeInverseIP,
   fprintf('Max error in inverse = %g\n', ...
           norm(S.BI*S.B - speye(size(G.cells.faces,1)), inf));
end

%--------------------------------------------------------------------------
%- Define meta data on pressure system component matrices -----------------
%
S.ip      = ipname;
S.type    = systemType;
S.intVert = opt.intVert;


% Inverse of innerproduct - functions

function W = ip_simple(a, v, K, C, N, computeInverseIP)
% These are not mutually inverse for general cells
t = 6* sum(diag(K))/size(K,2);
if computeInverseIP,
   Q  = orth(bsxfun(@times, C, a));
   U  = eye(length(a)) - Q*Q';
   d  = diag(a);
   W  = (N*K*N' + t*(d * U * d)) ./ v;
else
   Q  = orth(bsxfun(@rdivide, N, a));
   U  = eye(length(a)) - Q*Q';
   di = diag(1 ./ a);
   W  = (C * (K \ C'))./v + (v / t)*(di * U * di);
end

%--------------------------------------------------------------------------

function W = ip_tpf(a, v, K, C, N, computeInverseIP) %#ok
td = sum(C .* (N*K), 2) ./ sum(C.*C, 2);
%c2=sum(C.*C, 2);
%td = (sum(C .* (K*C'), 2)./c2).*(sum(C .* N, 2) ./ c2);
if computeInverseIP,
   W = diag(abs(td));
else
   W = diag(1 ./ abs(td));
end

%--------------------------------------------------------------------------

function W = ip_rt(a, v, K, C, N, computeInverseIP)
% Only for orthogonal normals.
%{
if size(K,1)==2,
   ix = [1,3,2,4];
   C  = C(ix,:);
   a  = a(ix,:);
end
%}
Ki = inv(K);
%Q  = kron(eye(size(K,2)), [1, 1]) ./ sqrt(2);
u  = diag(diag(Ki)) .* (v / 6);
di = diag(1 ./ a);
Q  = abs(N')*di/sqrt(2);
n = abs(N*N');
mask = n>eps*max(n(:));
W  = mask.*(C*(K\C')./v + di*Q'*u*Q*di);

if computeInverseIP, W = inv(W); end

%--------------------------------------------------------------------------

function W = ip_quasirt(a, v, K, C, N, computeInverseIP)
% Works for general grid, and coincide with ip_rt for orthogonal grids.
t = 6;
if computeInverseIP,
   Q  = orth(C);
   U  = eye(length(a)) - Q*Q';
   d  =  diag(diag(N*K*N'));
   W  = (N*K*N' + t*(U * d * U)) ./ v;
else
   Q  = orth(N);
   U  = eye(length(a)) - Q*Q';
   di = diag(1./(diag(N*K*N')));
   W  = (C * (K \ C'))./v + (v / t)*(U * di * U);
end

%--------------------------------------------------------------------------

function W = ip_dimfamily(a, v, K, C, N, computeInverseIP,t)
% Works for general grid, and coincide with ip_rt for orthogonal grids.
if computeInverseIP,
   Q  = orth(C);
   U  = eye(length(a)) - Q*Q';
   d  =  diag(diag(N*K*N'));
   W  = (N*K*N' + t*(U * d * U)) ./ v;
else
   Q  = orth(N);
   U  = eye(length(a)) - Q*Q';
   di = diag(1./(diag(N*K*N')));
   W  = (C * (K \ C'))./v + (v / t)*(U * di * U);
end
%--------------------------------------------------------------------------

function W = ip_quasitpf(a, v, K, C, N, computeInverseIP)
Q  = orth(C);
U  = eye(length(a)) - Q*Q';
W1 = N * K * N';

% ?? 2D ??
W  = (W1 + 2*U*diag(diag(W1))*U) ./ v;

if ~computeInverseIP, W = inv(W); end


%--------------------------------------------------------------------------
% Other private utility functions below.
%--------------------------------------------------------------------------

function s = id(s)
s = ['computeMimeticIPVE:', s];

%--------------------------------------------------------------------------

function [Type, ipfun, ipName, IP, invIP] = setup(opt)

ipName  = opt.InnerProduct;
Type    = opt.Type;

switch lower(Type),
   case 'hybrid',         [invIP, IP] = deal(true , false);
   case {'mixed','tpfa'}, [invIP, IP] = deal(false, true );
   case 'comp_hybrid',    [invIP, IP] = deal(true , true );
   otherwise,
      error(id('SystemType:Unknown'), ...
            'Unkown system type ''%s''.', Type);
end

ipNames = {'ip_simple', 'ip_tpf', 'ip_quasitpf', 'ip_rt', 'ip_quasirt','ip_dimfamily'};
if isempty(find(strcmp(ipName, ipNames), 1)),
   error(id('InnerProduct:Unknown'), ...
         'Unknown inner product ''%s''.', ipName);
end

if opt.verbose && strcmp(ipName, 'ip_rt'),
   disp(['NOTE: The mixed MFEM inner product is only valid for', ...
         ' Cartesian grids as given by cartGrid.'])
end

dispif(opt.verbose, 'Using inner product: ''%s''.\n', ipName);
switch lower(ipName),
   case 'ip_simple'  , ipfun = @ip_simple;
   case 'ip_tpf'     , ipfun = @ip_tpf;
   case 'ip_quasitpf', ipfun = @ip_quasitpf;
   case 'ip_rt'      , ipfun = @ip_rt;
   case 'ip_quasirt' , ipfun = @ip_quasirt;
   case 'ip_dimfamily' , ipfun = @(a, v, K, C, N, computeInverseIP) ...
         ip_dimfamily(a, v, K, C, N, computeInverseIP ,opt.inner_param);
   otherwise,
      % Can't happen.  'ipName' is already verified.
      error(id('InnerProduct:Unknown'), ...
            'Unknown inner product ''%s''.', ipName);
end




