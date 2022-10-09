%% Example script of diffusive time of flight
% Example of diffusive time of flight for a horizontal layer of the SPE10 
% model. The computeDTOF function is an implementation of the Fast Marching 
% Method (FMM) for diffusive time of flight as presented in 
% Zhang et al. 2013 (SPE 163637) 

%% User defined variables

%Fluid
viscosity           = 1*centi*poise;
density             = 1014*kilogram/meter^3; 
compressibility     = 4.4e-10;

%% Load necessary MRST modules
%%
mrstModule add spe10 heterogeneity incomp

%% load SPE10 - data, will attempt to download if not available locally
[G, W, rock] = getSPE10setup(30); % Loads layer
rock.perm(rock.poro==0,:) = 0; % To fix positive perm in deactivated cells.

%% Line-drive bottom to top
[nx, ny] = deal(G.cartDims(1), G.cartDims(2)); 
bottomCells = (1:nx)';
topCells    = (1:nx)' + nx*(ny-1);

%% Fluid structure
fluid       = initSingleFluid('mu' ,    viscosity, 'rho', density);
%% Calculating cell sizes
% Rectangular cuboids are assumed. In this case all cells are equal, and we
% could simply set cellSize = [6.096 3.048 0.6096] and get the same result.  
cellNo      = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
cf          = accumarray([cellNo, G.cells.faces(:,2)], G.cells.faces(:,1));
gc          = G.faces.centroids;
cellSize    = [gc(cf(:,2),1)-gc(cf(:,1),1) gc(cf(:,4),2)-gc(cf(:,3),2) ...
                gc(cf(:,6),3)-gc(cf(:,5),3)];
clear cf gc cellNo

%% Calculate diffusive time of flight
dtof=computeDTOF(G,rock,fluid,compressibility,bottomCells,cellSize);

%% Plotting diffusive time of flight
% The square is plotted as the units of "diffusive time of flight" as is
% square root of time.

figure(1), plotCellData(G, dtof.^2,...
 'Edgecolor', 'none', 'FaceAlpha', .6) ;
c=colorbar;
c.Label.String = 'Square of Diffusive Time of Flight';
axis off tight equal
