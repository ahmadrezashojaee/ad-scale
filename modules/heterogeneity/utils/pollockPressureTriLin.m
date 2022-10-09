function varargout = pollockPressureTriLin(G, state, rock, varargin)
% Trace streamlines in logically Cartesian grid using Pollock approximation.
% In addition to the regular Pollock approximation, pressures are
% supported by the use of trilinear interpolation (not recommended to use, 
% use pollockMod instead). 
% 
%
% SYNOPSIS:
%   [S,T,C,P] = pollockMod(G, state, rock)
%   [S,T,C,P] = pollockMod(G, state, rock, startpos)
%   [S,T,C,P] = pollockMod(G, state, rock, 'pn', pv, ...)
%   [S,T,C,P] = pollockMod(G, state, rock, startpos, 'pn', pv, ...)
%
% PARAMETERS:
%
%   G         - Cartesian or logically Cartesian grid.
%
%   state     - State structure with field 'flux'.
%
%   rock      - Rock structure with the field 'poro'.
%
% OPTIONAL PARAMETERS
%
%   positions - Matrix of size (N, 1) or (N, d+1), where d is the dimension
%               of the grid, used to indicate where the streamlines should
%               start.
%
%               If the size is (N, 1), positions contains the cell indices
%               in which streamlines should start. Each streamline is
%               started in the the local coordinate (0.5, 0.5, ...). To be
%               precise, this is the mean of the corner points, not the
%               centroid of the cell.
%
%               If the size is (N, d+1), the first column contains cell
%               indices, and the d next columns contain the local
%               coordinates at which to start streamlines.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%
%   substeps   - Number of substeps in each cell, to improve visual quality.
%               Default 5.
%
%   maxsteps   - Maximal number of points in a streamline.
%               Default 1000.
%
%   reverse    - Reverse velocity field before tracing.
%               Default false.
%
% RETURNS:
%
%  S      - Cell array of individual streamlines suitable for calls like
%           streamline(pollock(...)) and streamtube(pollock(...)).
%
%  T      - Time-of-flight of coordinate.
%
%  C      - Cell number of streamline segment, i.e, line segment between
%           two streamline coordinates.
%
%  P      - Interpolated pressure in each coordinate

% EXAMPLE:
%
%   [S,T,C,P] = pollockPressureTriLin(G, state, rock,startpos);
%
%  
% SEE ALSO: pollockMod()
%
%
%{
ORIGINAL COPYRIGHT FROM MRST: 

Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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

% Written by Jostein R. Natvig, SINTEF ICT, 2010.
%
% Modified by Asgeir Nyvoll, MSc student NTNU, 2018

   d = size(G.nodes.coords, 2);
   if mod(length(varargin),2)==0
      positions   = [(1:G.cells.num)', repmat(0.5, [G.cells.num, d])];
   else
      positions = varargin{1};
      if size(positions, 2) ==1
         positions = [positions, repmat(0.5, [size(positions, 1), d])];
      elseif size(positions, 2) ~= 1 + d
         error('Expected array of local positions of width 1 or 1+d.');
      end
      varargin  = varargin(2:end);
   end
   opt = struct('substeps', 5, 'maxsteps', 1000, 'reverse', false);
   opt = merge_options(opt, varargin{:});

   if opt.reverse
      state.flux = -state.flux;
   end

   [varargout{1:nargout}] = trace(G, state, positions, opt, rock);

end



% ========================================================================
function varargout = trace(G, state, pos, opt, rock)
   d              = size(G.nodes.coords, 2);
   numStreamlines = size(pos,1);
   assert(size(pos, 2) == d+1);

   if ~isfield(G, 'cellNodes')
      cn  = cellNodes(G);
      G.cellNodes = accumarray(cn(:,1:2), cn(:,3));
   end

   % Make array face fluxes for each cell in grid (Not outer).
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   cf     = G.cells.faces;
   pv     = poreVolume(G,rock);
   flux   = accumarray([cellNo, cf(:,2)], state.flux(cf(:,1)));
   velo   = flux./pv;
   velo(pv==0,:)=0;
   clear cf cellNo

   neighbors  = findNeighbors(G);

   nodePressure=nodePressures(G,state,rock);
   nodePressure(G.nodes.coords(:,2)==min(G.nodes.coords(:,2)))=5e7;
   nodePressure(G.nodes.coords(:,2)==max(G.nodes.coords(:,2)))=1e7;
   state.cellNodePressure=zeros(size(G.cellNodes));
   state.cellNodePressure(:,:)=nodePressure(G.cellNodes);
   clear nodePressure
   
   
   
   magic  = 1000;
   XYZ    = nan(numStreamlines, d, magic);
   T      = nan(numStreamlines, magic);
   C      = nan(numStreamlines, magic);
   P      = nan(numStreamlines, magic);
   active = true(numStreamlines, 1);


   % Store crossing coordinates of active streamlines
   [XYZ(active,:,1)] = globalCoordinate(G, pos(active,1), pos(active, 2:end));
   T(active, 1) = zeros(sum(active), 1);
   C(active, 1) = pos(active,1);
   P(active, 1) = localPressure(state,    pos(active,1), pos(active, 2:end));
   
   i = 2;
   while any(active)
      % Realloc
      if i+opt.substeps+1 > size(XYZ, 3)
         magic = max(magic, opt.substeps+1);
         XYZ   = cat(3, XYZ, nan(numStreamlines, d, magic));
         T     = cat(2, T,   nan(numStreamlines, magic));
         C     = cat(2, C,   nan(numStreamlines, magic));
         P     = cat(2, P,   nan(numStreamlines, magic));
      end
      current_cell = pos(active,1);

      % Take another pollock step
      [pos(active, :), t, xyz] = step(pos(active,:), velo, neighbors, opt.substeps);

      % Store crossing coordinates and, optionally, coordinates along curve
      % trajectory in cell of active streamlines
      for k=1:opt.substeps
         [XYZ(active, :, i+k-1)] = globalCoordinate(G, current_cell, xyz(:,:,k));
         P(active, i+k-1) = localPressure(state, current_cell, xyz(:,:,k));
      end
      T(active, i-1+(1:opt.substeps)) = repmat(t/opt.substeps, [1, opt.substeps]);
      C(active, i-1+(1:opt.substeps)) = repmat(pos(active, 1), [1, opt.substeps]);

      % Update active flag
      active(active)   =  pos(active,1) ~= current_cell;

      i = i+opt.substeps;
      if i > opt.maxsteps, break;end
   end

   %% Pack coordinates in list with streamlines separated by NaN.
   p = reshape(permute(XYZ, [3,1,2]), [], d);

   i = ~isnan(p(:,1));
   j = i|[true;i(1:end-1)];
   p = p(j,:);

   % Pack streamline coordinates in a cell array suitable for use with
   % Matlab streamline, i.e., as in 'streamline(pollock(G, resSol));'
   flag = isnan(p(:,1));
   ix = find(flag);
   dd  = diff([0;ix])-1;
   varargout{1} = mat2cell(p(~flag,:), dd, d);
   if nargout > 1
      T = reshape(T', [], 1);
      T = T(j);
      varargout{2} = mat2cell(T(~flag), dd, 1);
   end
   if nargout > 2
      C = reshape(C', [], 1);
      C = C(j);
      varargout{3} = mat2cell(C(~flag), dd, 1);
   end
   if nargout > 3
      P = reshape(P', [], 1);
      P = P(j);
      varargout{4} = mat2cell(P(~flag), dd, 1);
   end
end



% ========================================================================
function xyz = globalCoordinate(G, c, p)
% Compute global coordinate corresponding to local coorinate p in cells c
% p  - local positions == [xi,eta,zeta] in 3D
% c  -
%
   if numel(c)==1, p = reshape(p, 1, []); end
   % Compute node weight for quadrilateral or hexahedron
   d = size(G.nodes.coords, 2);
   w = ones(size(p,1), 2^d);
   for i=1:d
      mask        = logical(bitget((0:2^d-1)', i));
      w(:, mask)  = w(:, mask).* repmat(  p(:,i), [1, sum( mask)]);
      w(:,~mask)  = w(:,~mask).* repmat(1-p(:,i), [1, sum(~mask)]);
   end

   % Compute weighted average of corner points
   xyz = zeros(size(p,1), d);
   for i=1:d
      xi       = G.nodes.coords(:,i);
      xyz(:,i) = sum( w.*reshape(xi(G.cellNodes(c, :))', 2^d, [])', 2);
   end
end



%% ========================================================================
function [pos, tof, xyz] = step(pos, flux, neighbors, nsubsteps)
% Update pos array by computing new local coordinate and new cell.
% In addition, compute curve within cell.
%
%
%
%
   f = flux(pos(:,1),:);
   n = neighbors(pos(:,1),:);

   dims = size(pos, 2)-1;
   T    = nan(size(pos,1),dims);
   for i=1:dims
      T(:,i) = computeTime(pos(:,1+i), f(:,2*i-1:2*i));
   end
   [tof, dir] = min(T, [], 2);

   xyz = zeros(size(pos,1), dims, nsubsteps);
   d   = zeros(size(pos, 1), 1);
   for s=1:nsubsteps
      for i=1:dims
         t = tof*s/nsubsteps;
         [xyz(:,i,s), d(:,i)] = computePosition(pos(:,1+i), f(:,2*i-1:2*i), t);
      end
   end

   pos (:,2:end) = xyz(:,:,s);

   % Find direction to look up neighbor cell
   k  = 2*(dir-1)+d(sub2ind([numel(dir), 3], (1:numel(dir))', dir));
   t  = sub2ind(size(n), (1:numel(k))', k);

   % Update cell number if NOT at boundary.
   % IF at boundary, mark dir with NaN to avoid changing local coordinate
   % below.
   ind         = n(t)==0;
   pos(~ind,1) = n(t(~ind));
   dir (ind)   = nan;

   % Change local coordinate when moving to new cell
   k = sub2ind(size(d), (1:size(dir,1))', dir);
   k = k(~isnan(k));
   pos(numel(dir) + k ) = 2-d(k);
end


% ========================================================================
function N = findNeighbors(G)
% Build (n x 2*d) -array of neighbors for each cell in (Cartesian) grid G.
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   col    = 1 +   (cellNo == G.faces.neighbors(G.cells.faces(:,1), 1));
   c      = G.faces.neighbors(double(G.cells.faces(:,1)) + G.faces.num* (col-1));
   N      = accumarray([cellNo, G.cells.faces(:,2)], c);
end






% ========================================================================
function t = computeTime(xi, v)
% Compute time needed to reach xi=0 or xi=1 given velocities v=[v1,v2] at
% xi=0 and xi=1.  The formula is
%
%   t = xi/ui  or t = (1-xi)/ui,    if v1 = v2 = ui, and
%
%   t = 1/(v2-v1)*log(ue/ui),       otherwise
%
% where ui=v2*xi+v1*(1-xi) is the velocity at xi, and ue=v2 if ui>0 or
% ue=v1 if ui<0.
   tolerance = 100*eps;

   ui         = v(:,1) + xi.*diff(v, 1, 2);%(:,2)-v(:,1));
   ue         = v(:,    2);
   ue (ui<0)  = v(ui<0, 1);
   arg        = ue./ui;
   t          = inf(size(xi));

   % Linear velocity
   ind        = abs(diff(v, 1, 2)) > tolerance*abs(v(:,1));
   t(ind,:)   = 1./diff(v(ind,:), 1, 2).*log(arg(ind,:));

   % Constant velocity
   ds         = -xi;
   ds(ui > 0) = 1-xi(ui>0);
   t(~ind)    = ds(~ind)./ui(~ind);

   % nan happens for ui=ui=0
   t(arg<0 | isnan(arg))   = inf;
end


% ========================================================================
function [x, i] = computePosition(xi, v, t)
% Compute position at time t given start point xi and velocities v=[v1,v2].
%
%   x = xi + v*t,    if v is constant or
%
%   x = xi + (ui*exp((v2-v1)*t) - ui)/(v2-v1), otherwise
%
   tolerance = 100*eps;

   du        = diff(v, 1, 2);
   ui        = v(:,1) + xi.*du;
   i         = 1 + ~(ui<0);
   x         = inf(size(xi));

   ind       = abs(du) > tolerance*abs(v(:,1));

   % linear velocity
   x(ind)    = xi(ind) + ( ui(ind).*exp(du(ind).*t(ind)) - ui(ind))./du(ind);

   % Constant velocity
   x(~ind,:) = xi(~ind,:) + v(~ind, 1).*t(~ind, :);
   x(~ind & t==inf) = xi(~ind & t==inf);
end

function nodeP=nodePressures(G,state,rock)
%Calculates pressure in each node (grid cell corners)
cn=cellNodes(G);
nodeCells=accumarray([cn(:,3), cn(:,2)], cn(:,1));
weights=nodeCells;
cellPressures=weights;
active=nodeCells>0;
wx=weights; wy=weights; wz=weights;
wx(active)=G.cells.centroids(nodeCells(active),1);
wx=wx-G.nodes.coords(:,1);
wx(nodeCells==0)=0;

wy(active)=G.cells.centroids(nodeCells(active),2);
wy=wy-G.nodes.coords(:,2);
wy(nodeCells==0)=0;

wz(active)=G.cells.centroids(nodeCells(active),3);
wz=wz-G.nodes.coords(:,3);
wz(nodeCells==0)=0;


dims=size(rock.perm);
% Length and permeability scaled weighting
if(dims(2)==1)
weights(weights>0)=(abs(wx(active)).*rock.perm(nodeCells(active),1)...
    +abs(wy(active))+abs(wz(active))).*rock.perm(nodeCells(active))...
    ./(wx(active).^2+wy(active).^2+wz(active).^2);

else 
weights(weights>0)=(abs(wx(active)).*rock.perm(nodeCells(active),1)...
    +abs(wy(active)).*rock.perm(nodeCells(active),2)...
    +abs(wz(active)).*rock.perm(nodeCells(active),3))...
    ./(wx(active).^2+wy(active).^2+wz(active).^2);
    
end
weightsum=sum(weights,2);
cellPressures(active)=state.pressure(nodeCells(active));
nodeP=sum((weights.*cellPressures),2)./weightsum;

end

function streamlinePressure=localPressure(state,cell, pos)
%Trilinear interpolation
c000=state.cellNodePressure(cell,1);
c100=state.cellNodePressure(cell,2);
c010=state.cellNodePressure(cell,3);
c110=state.cellNodePressure(cell,4);
c001=state.cellNodePressure(cell,5);
c101=state.cellNodePressure(cell,6);
c011=state.cellNodePressure(cell,7);
c111=state.cellNodePressure(cell,8);

c00=c000.*(1-pos(:,1))+c100.*pos(:,1);
c01=c001.*(1-pos(:,1))+c101.*pos(:,1);
c10=c010.*(1-pos(:,1))+c110.*pos(:,1);
c11=c011.*(1-pos(:,1))+c111.*pos(:,1);

c0=c00.*(1-pos(:,2))+c10.*pos(:,2);
c1=c01.*(1-pos(:,2))+c11.*pos(:,2);

streamlinePressure=c0.*(1-pos(:,3))+c1.*pos(:,3);
end 




