%% Example script with Pollock algorithm on periodic grid 
% The pollockMod function is an expansion and modification of the regular
% pollock algorithm in MRST. This example is with a periodic grid in the
% x-direction of a single layer of the SPE10 model. 
%
% NB: Add 'Gp.faces.areas=Gp.faces.areas(1:Gp.faces.num);' at the end of 
% 'makePeriodicGridMulti3D' in the 'upscaling' module to fix problem with
% too many values in 'Gp.faces.areas'.

%% User-defined variables

%Streamlines
streamlinesPerCell  = 1;
nsubsteps           = 10;
linelength         = 'integral'; %Length estimation method

% Boundary Conditions 
pIn                 = 500*barsa();
pOut                = 100*barsa();

%Fluid
viscosity           = 1*centi*poise;
density             = 1014*kilogram/meter^3; 

% rateLimit: since the solver can give very small positive
% rates in cells that are enclosed by no flow cells, but still has positive
% permeability themselves. We want to avoid that streamlines can start in 
% such cells. (eps is the floating point relative accuracy of Matlab)

rateLimit           = 100*eps; % Cutoff for "active cells" 

%% Load necessary MRST modules

mrstModule add spe10 incomp heterogeneity upscaling
 
%% load SPE10 - data, will attempt to download if not available locally
[G, W, rock] = getSPE10setup(30); % Loads layer
is_pos = rock.poro > 0;
rock.poro(~is_pos) = min(min(rock.poro(is_pos)) / 100, 1.0e-6);

%% Periodic boundaries
% Periodic in x-direction
bcr{1}=pside([],G,'RIGHT',0);   bcl{1}=pside([],G,'LEFT',0);

%% Constant pressure boundaries: Linedrive bottom-to-top
% Constant pressure at y-direction boundaries (z-direction is no-flow)
[nx, ny] = deal(G.cartDims(1), G.cartDims(2)); 
bottomCells = (1:nx)';
topCells    = (1:nx)' + nx*(ny-1);
%Make periodic grid
[Gp, bcp]=makePeriodicGridMulti3d(G, bcl, bcr, {0});
bc = pside([],Gp,'YMin',pIn);
bc = pside(bc,Gp,'YMax',pOut);

 
%% Initialization of periodic grid
T       = computeTransGp(G, Gp, rock);
fluid   = initSingleFluid('mu' ,    viscosity, ...
                           'rho', density);
initState = initResSol(Gp, 0.0);
%% Solve incompTPFA to get steady-state pressure solution
 sol = incompTPFA(initState, Gp, T, fluid, 'bc', bc,'bcp',bcp);
    
%% Find flux in each cell to compare with the rate limit

sol.flux(isnan(sol.flux))=0;    %Fix NaN

cellNo = rldecode(1:Gp.cells.num, diff(Gp.cells.facePos), 2) .';
   cf     = Gp.cells.faces;
   flux= accumarray([cellNo, cf(:,2)], sol.flux(cf(:,1)));
   clear cf cellNo
   flux(:,2:2:end)=-1*flux(:,2:2:end);
cellFluxSum=sum(flux,2);
cellFlux=sum(flux.*(flux>0),2);


%% Compute start positions for streamlines
startCells = repmat(bottomCells,1,streamlinesPerCell);
startCells = reshape(startCells',[],1);
startPosLoc = [zeros(length(startCells),2), ones(length(startCells),1).*0.5];
startStep = 1/streamlinesPerCell; % Streamline spacing in start cell
                                  % in unit cell coordinates for pollock()  
first = startStep/2;       % Unit cell coordinate first streamline in cell                                    
last = 1-first;            % Unit cell coordinate first streamline in cell 
startPosLoc(:,1)=repmat((first:startStep:last)',length(bottomCells),1);
startPos = [startCells, startPosLoc]; 

% Finds active start cells
active=cellFlux(startCells)>rateLimit;

[S, T, C]=pollockMod(G, sol,rock,startPos(active,:), ...
    'substeps', nsubsteps ,'maxsteps', 4e5,'fluid',fluid, 'periodic', Gp, 'lineLength',linelength);

figure(1), plotCellData(G, log(convertTo(rock.perm(:,2),milli*darcy)),...
 'Edgecolor', 'none', 'FaceAlpha', .6) 
c=colorbar('Ticks',log([0.01,0.1,1,10,100,1000,10000]), 'TickLabels',...
{'0.01', '0.1','1','10','100','1000','10000'});
c.Label.String = 'Permeability [mD]';
title 'Streamlines on periodic grid'
axis off tight equal

% Plot streamlines if wanted (choose figure number from the plots above)
figure(1), h = streamline(S);
set(h, 'Color', 'black', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 3)
% increase z to plot streamlines above perm-field
for k =1:numel(h)
    h(k).ZData = h(k).ZData - 10;
end
