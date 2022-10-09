function [t, tC, bn, bnC] = getConnListAndBdyNodeWR2D(p, ny, na)
% Get connectivity list and boundary nodes of 2D well region (WR). The WR 
% is composed of a Cartesian region and two half-radial regions in xy 
% plane, which is used to connect the HW grid. For the Cartesian region, 
% the X axis extends along the well trajectory:
%     ---------> X
%    |    -----------------------------------
%  Y |    --------- Well trajectory ---------
%    V    -----------------------------------
%
% SYNOPSIS:
%   [t, tC, bn, bnC] = getConnListAndBdyNodeWR2D(p, ny, na)
%
% PARAMETERS:
%  p      - Points of the WR region, generated by 'pointsSingleWellNode'
%  ny     - The number of Cartesian cells in Y direction
%  na     - The number of angular cells in radial region
%
% RETURNS:
%  t      - Connectivity list of the whole well region
%  tC     - Connectivity list of the Cartesian region
%  bn     - Indices of outer boundary nodes of the whole well region
%  bnC    - Indices of outer boundary nodes of the Cartesian region
%
% SEE ALSO:
%  `VolumeOfInterest` `pointsSingleWellNode` `nearWellBoreModelingGrids`

    % Number of Cartesian nodes 
    nny  = ny + 1;
    nnx  = length(p);
    nnxy = nnx*nny;

    % Number of radial nodes 
    nr   = ny/2;
    nnra = nr*(na-1);

    % Cartesian node indices, reshaped to 2D
    ndC = reshape((1:nny*nnx)', nny, nnx);

    % Connectivity list for the radial nodes, Heel
    ndR1 = reshape(nnxy + (1:nnra)', nr,  na-1);
    ndR1 = [(1:nr)', ndR1, (nny:-1:nr+2)'];
    ndR1 = [ndR1; (nr+1) * ones(1, size(ndR1,2))];

    % Connectivity list for the radial nodes, Toe
    ndR2 = reshape(nnxy + nnra + (1:nnra)', nr,  na-1);
    d  = (nnx-1)*nny;
    ndR2 = [(1:nr)' + d, ndR2, (nny:-1:nr+2)' + d];
    ndR2 = [ndR2; (nr+1 + d) * ones(1, size(ndR2,2))];

    % Connectivity list for the Cartesian nodes
    tC  = makeConnListFromMat(ndC);
    tR1 = makeConnListFromMat(ndR1);
    tR2 = makeConnListFromMat(ndR2);

    % Combine the connectivity lists
    t = [tC; tR1; tR2];
    idxTri = cellfun(@(x)length(x) ~= length(unique(x)), t);
    t(idxTri) = cellfunUniOut(@unique, t(idxTri));

    % Get boundary nodes, with radial region
    bn = [ndC(1,:)'; ndR2(1, 2:end-1)'; ndC(end,end:-1:1)'; ...
        ndR1(1, end-1:-1:2)'];

    % Get boundary nodes, without radial region
    bnC = [ndC(1,:)'; ndC(2:end-1,end); ndC(end,end:-1:1)';...
        ndC(end-1:-1:2,1);];
end