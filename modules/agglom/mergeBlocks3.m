function p = mergeBlocks3(p, G, I, N, flist, varargin)
%Amalgamation 'MERGE' primitive adapted to fault information
%
% SYNOPSIS:
%   p = mergeBlocks3(p, G, I, N, flt)
%   p = mergeBlocks3(p, G, I, N, flt, 'pn1', pv1, ...)
%
% PARAMETERS:
%   p   - Partition vector.  Possibly created by function
%         'applySuccessivePart' or some other partitioning algorithm.  The
%         input partition vector should not contain any disconnected
%         blocks.  Function 'processPartition' will split such blocks.
%
%   G   - Grid structure.
%
%   I   - Indicator function structure.  Each indicator function,
%         represented as a non-negative scalar for each active cell in the
%         grid 'G', is interpreted as a density.  The following indicators
%         must be present as individual structure fields:
%
%            - volume -- Block volume indicator (e.g., porosity)
%            - flow   -- Block flow indicator (e.g., time-of-flight)
%
%         The algorithm will generally merge a candidate block to the its
%         (feasible) neighbouring block that most closely matches its own
%         block flow.
%
%   N   - Algorithm controlling parameters.  Structure that features the
%         following fields
%
%            - volumeLow --
%                 Lower bound on total block volume.  The algorithm will
%                 merge blocks that violate the criterion
%
%                    I.volume(B) |B| >= (N.volumeLow / n) I.volume(G) |G|
%
%            - flowHigh --
%                 Upper bound on total block flow.  While merging blocks
%                 that violate the lower block volume criterion, it will
%                 attempt to uphold the upper bound flow criterion
%
%                    I.flow(B) |B| <= (N.flowHigh / n) I.flow(G) |G|
%
%            - blocksHigh --
%                 Upper bound on total number of fine-scale cells in a
%                 block.  When merging blocks that violate the lower block
%                 volume criterion, the algorithm will attempt to uphold
%                 the criterion
%
%                    #cells(B) <= N.blocksHigh
%
%                 OPTIONAL.  Treated as INF (inactive) if not present.
%
%   flt - M-element structure array that defines faults.  Each structure
%         element must define the following fields:
%
%            - hard -- Numeric array of faces (face indices) that
%                 constitute the fault's connections.
%
%            - soft -- Numeric array of faces (face indices) that extend
%                 the fault towards the boundary.
%
% OPTIONAL PARAMETERS (supplied in 'key'/value pairs ('pn'/pv ...)):
%   cfac - A relative factor at which the upper bound(s) are turned into
%          hard constraints.  Default: cfac = INF (do not turn criteria into
%          hard constraints).
%
%   merge_vert -
%          LOGICAL (Boolean) flag indicating whether or not to merge blocks
%          in the vertical direction.  Default value: merge_vert = TRUE (do
%          merge blocks in vertical direction).
%
% RETURNS:
%   p - Updated partition vector.  Typically contains fewer blocks than the
%       input partition vector.  None of the resulting blocks should
%       violate the criterion lower bound on block volumes, but some of the
%       blocks may violate the upper bounds on total block flow and/or
%       total number of cells in a block.
%
% SEE ALSO:
%   `mergeBlocks2`, `refineBlocks`, `processPartition`.

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

   opt = struct('cfac', inf, 'merge_vert', true);
   opt = merge_options(opt, varargin{:});

   bI    = block_indicators(G, p, I);
   bN    = define_neighbourship(G, p, flist, opt.merge_vert);
   bound = compute_bounds(G.cells.num, bI, N);

   mrg = mergeBlocksCore(bN, bI, bound, opt.cfac);

   p = compressPartition(mrg(p));
   q = processPartition(G, p);
   if ~all(p == q),
      warning(msgid('Part:Disconnected'), ...
             ['Resulting partition is disconnected. ', ...
              'Repairing before continuing.']);
      p = q;
   end

   if any(accumarray(p, I.volume .* G.cells.volumes) < bound.volumeLow),
      warning(msgid('LBnd:Violated'), ...
             ['Some blocks still violate lower (volume) ', ...
              'bound after merging.']);
   end
end

%--------------------------------------------------------------------------

function bI = block_indicators(G, p, I)
   accum = sparse(p, 1 : G.cells.num, 1) * ...
           [ (I.volume .* G.cells.volumes), ...
             (I.flow   .* G.cells.volumes), ...
             ones([G.cells.num, 1]) ];

   bI = struct('Vol', accum(:,1), ...
               'Flw', accum(:,2) ./ accum(:,1), ...
               'Num', accum(:,3));
end

%--------------------------------------------------------------------------

function bN = define_neighbourship(G, p, flist, vertical)
   neigh_kind = 'Topological';  % Include explicit NNCs

   % Include boundary connections.  Needed for consistent indexing using
   % 'flist.hard'
   incBdry = true;

   N = getNeighbourship(G, neigh_kind, incBdry);

   pick = true([size(N, 1), 1]);

   if ~vertical,
      vert_conn_p      = false([max(G.cells.faces(:,2)), 1]);
      vert_conn_p(5:6) = true;

      is_vert = accumarray(G.cells.faces(:,1), ...
                           vert_conn_p(G.cells.faces(:,2)), ...
                           [size(N,1), 1], @any);

      pick = pick & ~is_vert;  clear is_vert vert_p
   end

   if ~isempty(flist),
      pick = pick & eliminate_fault_conn(N, p, flist);
   end

   bN = blockNeighbourship(N(pick, :), p);
end

%--------------------------------------------------------------------------

function bound = compute_bounds(nc, bI, N)
   bound = struct( ...
      'volumeLow', N.volumeLow  * sum(          bI.Vol) / nc, ...
      'flowHigh' , N.flowHigh   * sum(bI.Flw .* bI.Vol) / nc, ...
      'cellHigh' , inf);

   if isfield(N, 'blockHigh'),
      bound.cellHigh = N.blockHigh;
   end
end

%--------------------------------------------------------------------------

function pick = eliminate_fault_conn(N, p, flist)
   % Eliminate hard fault connections as well as (soft) extensions of those
   % hard faults that happen to affect the blocks intersected by the hard
   % faults.  Extensions in blocks unaffected by hard faults are treated as
   % regular connections.

   p1 = [ 0 ; reshape(p, [], 1) ];
   nh = cellfun('prodofsize', { flist.hard });  % #hard fault conns
   ns = cellfun('prodofsize', { flist.soft });  % #soft fault conns

   % hard = [ block , (hard) fault ID ]
   %
   % Note "2 * nh" to account for block pairs.
   hard = [ p1(reshape(N(vertcat(flist.hard), :) .', [], 1) + 1), ...
            rldecode(1 : numel(flist), 2 * nh, 2) .' ];

   hard = unique(hard, 'rows'); % Expensive.  Hope 'hard' isn't too large.

   % soft = concatenation of soft extension faces, one copy for each row in
   % 'hard'.  Potentially large.
   soft = vertcat(flist(hard(:, 2)).soft);

   % i = 'soft' indices for which one of the connecting blocks is affected
   % by the corresponding hard fault.
   i = any(bsxfun(@eq, p1(N(soft, :) + 1), ...
                  rldecode(hard(:, 1), ns(hard(:, 2)))), 2);

   % Mark hard faults and block-restricted soft extensions for exclusion.
   pick = true([size(N, 1), 1]);
   pick([soft(i) ; vertcat(flist.hard)]) = false;
end
