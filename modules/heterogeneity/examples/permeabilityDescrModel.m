%% Example script calculating permeability descriptors
% Model described in Nyvoll (2017) and Nyvoll (2018), which is a 
% specialization project and a master thesis, respectively. (Both currently 
% unpublished as of June 2018) 

%% User-defined variables
%Streamlines
streamlinesPerCell  = 20;
nsubsteps           = 10;
linelength         = 'straight';
% Boundary Conditions
pIn                 = 500*barsa();
pOut                = 100*barsa();

%Fluid
viscosity           = 1*centi*poise;
density             = 1014*kilogram/meter^3; 

% rateLimit: since the solver can give very small positive
% rates in cells that are enclosed by no flow cells, but still has positive
% permeability themselves. eps is the floating point relative accuracy
% of Matlab

rateLimit           = 100*eps; % Cutoff for "active cells" 


%% Load necessary MRST modules

mrstModule add spe10 incomp heterogeneity


%% load SPE10 - data, will attempt to download if not available locally
%%
layerNum = 1;  % choose layer 
[G, W, rock] = getSPE10setup(layerNum); % Loads layer
rock.perm(rock.poro==0,:)=0; % To fix positive perm in deactivated cells.
 
%% Constant pressure boundaries: Linedrive bottom-to-top
%%
[nx, ny] = deal(G.cartDims(1), G.cartDims(2)); 
bottomCells = (1:nx)';
topCells    = (1:nx)' + nx*(ny-1);
bc = pside([],G,'YMin',pIn);
bc = pside(bc,G,'YMax',pOut);
 
%% Computing Transmissibilities, fluid and create initial state
%%
T       = computeTrans(G, rock);
fluid   = initSingleFluid('mu' ,    viscosity, ...
                           'rho', density);
 initState = initResSol(G, 0.0);
%% Solve incompTPFA to get steady-state pressure solution
%%

 sol = incompTPFA(initState, G, T, fluid, 'bc', bc);
    
%% Find cell fluxes to compare with rate limit
sol.flux(isnan(sol.flux))=0;    %Fix NaN
cellNo = rldecode(1:G.cells.num, diff(G.cells.facePos), 2) .';
   cf     = G.cells.faces;
   flux= accumarray([cellNo, cf(:,2)], sol.flux(cf(:,1)));
   clear cf cellNo
   flux(:,2:2:end)=-1*flux(:,2:2:end);
cellFluxSum=sum(flux,2);
cellFlux=sum(flux.*(flux>0),2);


%% Compute start positions for streamlines
% Using equal spacing at each inlet face.
startCells = repmat(bottomCells,1,streamlinesPerCell);
startCells = reshape(startCells',[],1);
startPosLoc = [zeros(length(startCells),2), ones(length(startCells),1).*0.5];
startStep = 1/streamlinesPerCell; % Streamline spacing in start cell
                                  % in unit cell coordinates for pollockMod()  
first = startStep/2;       % Unit cell coordinate first streamline in cell                                    
last = 1-first;            % Unit cell coordinate first streamline in cell 
startPosLoc(:,1)=repmat((first:startStep:last)',length(bottomCells),1);
startPos = [startCells, startPosLoc]; 

% Finds active start cells
active=cellFlux(startCells)>rateLimit;

%% Modified pollock method
% Uses all possible outputs of pollockMod, including (in order): 
% coordinates, times of flight, current cell, streamline step lengths,
% velocity vectors and pressure gradients based on 1D darcy flow along each
% streamline.
[S, T, C,LS,VXYZ,GradP]=pollockMod(G, sol,rock,startPos(active,:), ...
    'substeps', nsubsteps ,'maxsteps', 4e5,'fluid',fluid, 'lineLength',linelength);

%% Various calculations
% Rate for each streamtube around streamlines:
Qs=ones(size(S)).*flux(startCells(active),3)./streamlinesPerCell;

% Total rate
Q=sum(Qs);

%Bulk volume
V=sum(G.cells.volumes);

% Effective bulk volume
Ve=sum(G.cells.volumes(cellFlux>rateLimit));

% Pore Volume
OmegaS=sum(sum(poreVolume(G,rock)));

% Find number of streamlines started in active cells:
nstreamlines=length(S); 
BS=zeros(nstreamlines,1);
CS=BS;
TauS2=BS;
% Model length in main flow direction (y-direction). 
Ltot=max(G.faces.centroids(:,2))-min(G.faces.centroids(:,2));

%% Streamline/streamtube permeability descriptors
for i=1:nstreamlines
   [BS(i), CS(i), TauS2(i)]=streamtubePermDesc(GradP{i}, LS{i}, fluid, 670.56);
   % BS streamline conductance, CS streamline constriction factor and TauS2
   % streamline tortuoisty factor.
end
tof=zeros(length(T),1);
for i = 1 : length(T)
    % Time of flight /residence time of each streamline
    tof(i,1)=sum(T{i}(:)); 
end

%% Effective layer permeability descriptors and variance
[Be,Ce,Te2,varTs2,varInvCs]= ...
    effPermDesc(BS,CS,TauS2, Qs, Ve);
 
%% Permeability
%ke calculated with Darcy over full layer
keff=-Q*fluid.properties(1)*670.56/...
    (-4e7*sum(G.faces.areas(13421:13480)));
%Average k_h in V (k_x=k_y=k_h)
kavg=mean(rock.perm(:,2));
%Average k_h in Ve (for spe10: Should equal Be)
kavgact=mean(rock.perm(cellFlux>rateLimit,2));
%ke calculated with effective descriptors from streamlines
keffs=Te2*Be*Ve/(Ce*V);
