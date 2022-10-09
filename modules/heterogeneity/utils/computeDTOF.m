function DTOF=computeDTOF(G,rock,fluid,comp,startCells,cellSize)
% Compute diffusive time of flight with Fast Marching Method (FMM). See
% Zhang et al. 2013 (SPE 163637). The diffusive time of flight is related
% to the propagation of a pressure front, similar to the depth of
% investigation.
% 
%
% SYNOPSIS:
%   DTOF = computeDTOF(G,rock,fluid,comp,startCells,cellSize)
%
% PARAMETERS:
%
%   G           - Cartesian or logically Cartesian grid.
%
%   rock        - Rock structure with the fields 'poro' and 'perm'.
%
%   fluid       - MRST fluid structure where fluid.properties(1)=viscosity
%
%   comp        - Total compressibility. Vector of size (1,1) or (N,1) 
%                 where N is the number of grid cells.
%
%   startCells  - Cells where DTOF=0;
%
%   cellSize    - Vector of size (1,d) or (N,d) where d is dimension of the
%                 grid and N is the number of cells.
%
% RETURNS:
%
%  DTOF   - Diffusive time of flight 
%
%  Written by Asgeir Nyvoll, MSc student NTNU, 2018

nbrs=findNeighbors(G);
DTOF=NaN(G.cells.num,1);
voxelState=zeros(G.cells.num,1); % 0 is far, 1 is evaluated, 2 is accepted
DTOF(startCells)=0;              % initial
voxelState(startCells)=2;        % accept initial
startNbrs=unique(nbrs(voxelState==2,:)); 
diffu=rock.perm./(rock.poro.*comp.*fluid.properties(1));
slowness=cellSize./sqrt(diffu)./2;
gdim=size(nbrs,2)/2;
if gdim==3
    if all(nbrs(:,5:6)==0)
        gdim=2; %Basically 2D
    end
end
%Preparing first set of neighbors
for n = startNbrs(:).'
    if n>0 && voxelState(n)==0
    voxelState(n)=1;
    DTOF(n)=solveEikonalDTOF(DTOF,voxelState,slowness,nbrs(n,:),n,gdim);
    end
end
% Accepting first neighbor
acceptedValue=min(DTOF(voxelState==1));      % Lowest DTOF is accepted
acceptedNode=find(DTOF==acceptedValue,1);    % Node of accepted DTOF
voxelState(acceptedNode)=2;                  % Update status
unacceptedNbrs=unique(nbrs(acceptedNode,:)); % Neighbors of accepted
unacceptedNbrs=unacceptedNbrs(unacceptedNbrs>0); % Remove boundary
% Remove accepted neighbors
unacceptedNbrs=unacceptedNbrs(voxelState(unacceptedNbrs)~=2);
voxelState(unacceptedNbrs)=1; %update state of unaccepted neighbors

% Iterate until all connected cells are reached
while ~isempty(acceptedNode) && any(voxelState~=2)
    for n=unacceptedNbrs(:).'
        % Update DTOF of unaccepted neighbors of latest accepted node
        DTOF(n)=solveEikonalDTOF(DTOF,voxelState,slowness,nbrs(n,:),n,gdim);
    end  
    % Accept next node, find new neighbors, update states
    acceptedValue=min(DTOF(voxelState==1));
    acceptedNode=find(DTOF==acceptedValue);
    voxelState(acceptedNode)=2;
    unacceptedNbrs=unique(nbrs(acceptedNode,:));
    unacceptedNbrs=unacceptedNbrs(unacceptedNbrs>0);
    unacceptedNbrs=unacceptedNbrs(voxelState(unacceptedNbrs)~=2);
    voxelState(unacceptedNbrs)=1;
end
end

function tempDTOF=solveEikonalDTOF(DTOF,voxelState,slowness,nbrs,n,gdim)
%Solves the Eikonal equation with the approximation from Zhang et al.
%
tempDTOF=NaN;

% Evaluate diffusive time of flight from each corner. Smallest value is
% valid
if gdim==3
    for i=1:2
        for j=3:4
            for k=5:6
                tempDTOF=min(singleTestSolver(...
                    DTOF,voxelState,slowness,nbrs,n,[i j k]),tempDTOF);
            end
        end
    end
elseif gdim==2
    for i=1:2
        for j=3:4
                tempDTOF=min(singleTestSolver(...
                DTOF,voxelState,slowness,nbrs,n,[i j]),tempDTOF);
        end
    end
else
    error('Neither 2D nor 3D')
    %Should "never" happen
end
end


function dtof=singleTestSolver(DTOF,voxelState,slowness,nbrs,n,testCase)
%Solves Equations (9) and (10) in Zhang et al. (2013) for one corner
coeff=[0 0 0]; % coeff = [a,b,c] in quadratic equation a*t^2+b*t+c
for i=1:length(testCase)
if (nbrs(testCase(i))>0) %Check that it is not a boundary
    direction=ceil(testCase(i)/2); 
    if isnan(slowness(n,direction)) 
        dtof=NaN;
        return
    end
    if (voxelState(nbrs(testCase(i)))==2 && ...
            ~isnan(slowness(nbrs(testCase(i)),direction)))
        coeff=coeff+1/(slowness(nbrs(testCase(i)),direction)+...
            slowness(n,direction))^2.*...
            [1, -2*DTOF(nbrs(testCase(i))),(DTOF(nbrs(testCase(i))))^2];
    end
end
end
if all(coeff==0)
    dtof=NaN; return
else
    coeff=coeff-[0 0 1];
    % Solve quadratic equation 
    dtof=((-coeff(2)+sqrt(coeff(2)^2-4*coeff(1)*coeff(3)))/(2*coeff(1)));
    return
end
end


function N = findNeighbors(G)
% Build (n x 2*d) -array of neighbors for each cell in (Cartesian) grid G.
   cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2)';
   col    = 1 +   (cellNo == G.faces.neighbors(G.cells.faces(:,1), 1));
   c      = G.faces.neighbors(double(G.cells.faces(:,1)) + G.faces.num* (col-1));
   N      = accumarray([cellNo, G.cells.faces(:,2)], c);
end

