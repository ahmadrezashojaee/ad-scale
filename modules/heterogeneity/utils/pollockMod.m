function varargout = pollockMod(G,state, rock, varargin)
% Trace streamlines in logically Cartesian grid using Pollock approximation.
% In addition to the regular Pollock approximation, pressure gradients are
% supported. 
% 
%
% SYNOPSIS:
%   [S,T,C,L,V,GP] = pollockMod(G, state, rock)
%   [S,T,C,L,V,GP] = pollockMod(G, state, rock, startpos)
%   [S,T,C,L,V,GP] = pollockMod(G, state, rock, 'pn', pv, ...)
%   [S,T,C,L,V,GP] = pollockMod(G, state, rock, startpos, 'pn', pv, ...)
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
%   periodic   - Grid with periodic boundary conditions.
%               Default not periodic boundary conditions.
%     
%   fluid      - MRST fluid structure where fluid.properties(1)=viscosity
%               Default viscosity: 1 cP
%
%   lineLength - Method for streamline length calculation. Options:
%                'straight' (straight line between plots) and 'integral'.
%                Default: 'straight'
%
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
%  L      - Streamline lengths
%
%  V      - Velocity vectors
%
%  GP     - Pressure gradients
%
% EXAMPLE:
%
%   [S,T,C,L,V,GP] = pollockMod(G, x, rock,startpos,'fluid',fluid);
%
%   streamline(S);
% SEE ALSO: pollock()
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
   
   opt = struct('substeps', 5, 'maxsteps', 1000, 'reverse', false, 'periodic',...
       [], 'fluid', [],'lineLength','straight');
   opt = merge_options(opt, varargin{:});
	
% Check if fluid properties are given for Darcy calculation of pressure gradient
    if isempty(opt.fluid) 
	opt.fluid.properties(1)=1*centi*poise; %default value
    end

    
% Check that streamline length is calculated with integral of inverse 
% velocity in case of periodic grid.
	if ~isempty(opt.periodic) && ~isequal(opt.lineLength,'integral')
	warning('Length calculation changed: Periodic grids only support integral') 
	opt.lineLength='integral';
	end
	
   if size(state.flux, 2) > 1
        state.flux = sum(state.flux, 2);
   end


   if opt.reverse
      state.flux = -state.flux;
   end

   [varargout{1:nargout}] = trace(G, state, positions, rock, opt);

end



% ========================================================================
function varargout = trace(G, state, pos, rock, opt)

if isempty(opt.periodic)
	%State grid and coordinate grid are equal when grid is not periodic.
	Gp=G;
else
	Gp=opt.periodic;
end

   d              = size(G.nodes.coords, 2);
   numStreamlines = size(pos,1);
   assert(size(pos, 2) == d+1);

   if ~isfield(G, 'cellNodes')
      cn  = cellNodes(G);
      G.cellNodes = accumarray(cn(:,1:2), cn(:,3));
   end
   

   % Make array face fluxes for each cell in grid
   cellNo = rldecode(1:Gp.cells.num, diff(Gp.cells.facePos), 2) .';
   cf     = Gp.cells.faces;
   pv     = poreVolume(Gp,rock);
   flux   = accumarray([cellNo, cf(:,2)], state.flux(cf(:,1)));
   
   %Particle velocity in unit cell
   unitVelo   = flux./pv;
   %Remove NaN cases:
   unitVelo(pv==0,:)=0;
   
   %Actual velocity
   velo1=state.flux./Gp.faces.areas;
   velo   = accumarray([cellNo, cf(:,2)],velo1(cf(:,1)))./rock.poro;
   velo(pv==0,:)=0;
   clear cf cellNo pv flux velo1
  
   neighbors  = findNeighbors(Gp);

   magic  = 1000;
   XYZ    = nan(numStreamlines, d, magic);
   VXYZ   = nan(numStreamlines, d, magic);
   T      = nan(numStreamlines, magic);
   C      = nan(numStreamlines, magic);
   Ls     = nan(numStreamlines, magic);
   GradP  = nan(numStreamlines, magic);
   
   active = true(numStreamlines, 1);

   % Store initial values
   [XYZ(active,:,1)] = globalCoordinate(G, pos(active,1), pos(active, 2:end));
   [VXYZ(active,:,1)] = globalVelocity(pos(active,1), pos(active, 2:end), velo);
   T(active, 1) = zeros(sum(active), 1);
   Ls(active, 1) = zeros(sum(active), 1);
   C(active, 1) = pos(active, 1);
   GradP(active, 1) = zeros(sum(active), 1); 

   i = 2;
   while any(active)
      % Realloc
   
      if i+opt.substeps+1 > size(XYZ, 3)
         magic = max(magic, opt.substeps+1);
         XYZ   = cat(3, XYZ,    nan(numStreamlines, d, magic));
         VXYZ  = cat(3, VXYZ,   nan(numStreamlines, d, magic));
         T     = cat(2, T,      nan(numStreamlines, magic));
         C     = cat(2, C,      nan(numStreamlines, magic));
         Ls    = cat(2, Ls,     nan(numStreamlines, magic));
         GradP = cat(2, GradP,  nan(numStreamlines, magic));  
      end
      current_cell = pos(active,1);
      C(active, i-1+(1:opt.substeps)) = repmat(pos(active, 1), [1, opt.substeps]);
      % Take another pollock step
      [pos(active,:), t, xyz, VXYZ(active, :, i:i+opt.substeps-1), ...
           Ls(active,i:i+opt.substeps-1), ...
           GradP(active,i:i+opt.substeps-1)] = ...
           step(rock, opt.fluid, pos(active,:), unitVelo, ...
           velo, neighbors, opt.substeps);
      T(active, i-1+(1:opt.substeps)) = repmat(t/opt.substeps, [1, opt.substeps]);

      %Store coordinates
      for k=1:opt.substeps
      
        %Global coordinate (using coordinate grid G for periodic grid Gp)
        [XYZ(active, :, i+k-1)] = globalCoordinate(G, current_cell, xyz(:,:,k));
        
        %Streamline step length (various options)
        if isequal(opt.lineLength,'straight')
            %Straight line between neighboring coordinates
            [Ls(active, i+k-1)] = ...
                (sqrt(sum((XYZ(active, :, i+k-1)-XYZ(active, :, i+k-2)).^2,2)));

        elseif isequal(opt.lineLength,'integral') 
           continue; %Already done in step method
        else
            error('Invalid streamline length method. Choose straight or integral');
        end
	
      end
      
     
      % Update active flag
      active(active)   =  pos(active,1) ~= current_cell;

      i = i+opt.substeps;
      %Break if reaching maxsteps. Display warning.
      if i > opt.maxsteps,warning('Maxsteps reached'), break;end
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
   % Pack times of flight.
   if nargout > 1
      T = reshape(T', [], 1);
      T = T(j);
      varargout{2} = mat2cell(T(~flag), dd, 1);
   end
   % Pack cells.
   if nargout > 2
      C = reshape(C', [], 1);
      C = C(j);
      varargout{3} = mat2cell(C(~flag), dd, 1);
   end
     
   % Pack step lengths.
   if nargout > 3
      Ls = reshape(Ls', [], 1);
      Ls = Ls(j);
      varargout{4} = mat2cell(Ls(~flag), dd, 1);
   end
   
   % Pack velocities.
   if nargout>4
   v = reshape(permute(VXYZ, [3,1,2]), [], d);

   i = ~isnan(v(:,1));
   j = i|[true;i(1:end-1)];
   v = v(j,:);

   flag = isnan(v(:,1));
   ix = find(flag);
   dd  = diff([0;ix])-1;
   varargout{5} = mat2cell(v(~flag,:), dd, d);
   end
   
   % Pack pressures gradients.
   if nargout > 5
     GradP = reshape(GradP', [], 1);
     GradP = GradP(j);
     varargout{6} = mat2cell(GradP(~flag), dd, 1);
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


function vxyz = globalVelocity(c, p, v)
% Return velocity vector at each substep.

d = size(v,2)/2;
w = size(c,1);
vxyz=zeros(w,d);

    for i=1:d
        %Linear interpolation internally in unit cell
        vxyz(:,i) = (v(c,2*i)-v(c,2*i-1)).*p(:,i)+v(c,2*i-1);

    end
end
%% ========================================================================
function [pos, tof, xyz, vxyz, ds, gp] = ...
    step(rock, fluid, pos, unitflux, flux, neighbors, nsubsteps)
% Update pos array by computing new local coordinate and new cell.
% In addition, compute curve within cell, step lengths, time of flight 
% and pressure gradients.
%
%
   uf = unitflux(pos(:,1),:);   % velocity unit cell
   f  = flux(pos(:,1),:);       % velocity 
   n = neighbors(pos(:,1),:);
   dims = size(pos, 2)-1;
   T    = nan(size(pos,1),dims);
   for i=1:dims
      T(:,i) = computeTime(pos(:,1+i), uf(:,2*i-1:2*i)); 
   end
   [tof, dir] = min(T, [], 2);
   % Compute positions, velocities, lengths and pressure gradients
   [xyz, d, vxyz, ds, gp] = ...
       computePosVelLenGrad(rock, fluid, pos, uf, tof, f, nsubsteps);
   pos (:,2:end) = xyz(:,:,end);

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
function [xyz, d,vxyz,ds,gp] = ...
    computePosVelLenGrad(rock, fluid, pos, uf, tof, f, nsubsteps)
% Compute position at time t given start point xi and velocities v=[v1,v2].
%
%   x = xi + v*t,    if v is constant or
%
%   x = xi + (ui*exp((v2-v1)*t) - ui)/(v2-v1), otherwise
%
   dims = size(pos, 2)-1;
   nel = size(pos, 1);
   xyz = zeros(nel, dims, nsubsteps);
   vxyz= xyz;
   ds  = zeros(nel, nsubsteps);
   gp  = ds;
   d   = zeros(nel, 1);
   du  = zeros(nel,dims); 
   dv  = du;
   ui  = du;
   vi  = du;
   dt= tof/nsubsteps;
    
   for i=1:dims
   du(:,i)= diff(uf(:,2*i-1:2*i), 1, 2); 
   dv(:,i)= diff(f(:,2*i-1:2*i), 1, 2);  
   ui(:,i)= uf(:,2*i-1) + pos(:,1+i).*du(:,i); %unit cell velocity
   vi(:,i)= f(:,2*i-1) + pos(:,1+i).*dv(:,i);  %velocity
   d(:,i)    = 1 + ~(ui(:,i)<0);
   end
  
for s=1:nsubsteps
    t = s*dt; %residence time in cell
      for i=1:dims
        tolerance = 100*eps;

        xyz(:,i,s)= inf(size(pos(:,1+i)));
         
        ind       = abs(du(:,i)) > tolerance*abs(uf(:,2*i-1)); 
        
        % linear velocity
        xyz(ind,i,s)   = pos(ind,1+i) + ...
            ( ui(ind,i).*exp(du(ind,i).*t(ind)) - ui(ind,i))./du(ind,i);
        vxyz(ind,i,s)  =(vi(ind,i)).*exp(du(ind,i).*t(ind));
        
        
        % Constant velocity
        xyz(~ind,i,s) = pos(~ind,1+i) + uf(~ind,2*i-1).*t(~ind, :);
        xyz(~ind & t==inf,i,s) = pos(~ind & t==inf,1+i);
        vxyz(~ind,i,s)=f(~ind,2*i-1);
                
      end
      
      if s==1
          vp=sqrt(sum(vi.^2,2));
          v=sqrt(sum(vxyz(:,:,s).^2,2)); %velocity at end of step
          %mean perm
          k=(sqrt(sum((rock.perm(pos(:,1),:).*vi./vp).^2,2))...
             +sqrt(sum((rock.perm(pos(:,1),:).*vxyz(:,:,s)./v).^2,2)))./2; 
      else
          vp=v;
          v=sqrt(sum(vxyz(:,:,s).^2,2)); %velocity at end of step
          %mean perm
          k=(sqrt(sum((rock.perm(pos(:,1),:).*vxyz(:,:,s-1)./vp).^2,2))...
             +sqrt(sum((rock.perm(pos(:,1),:).*vxyz(:,:,s)./v).^2,2)))./2; 
      end
         
          %curve length
          ds(:,s)=dt.*(vp+v)./2;
         
          %mean pressure gradient
          gp(:,s)= -((v+vp)./2.*rock.poro(pos(:,1))./k).*fluid.properties(1);
end
end








