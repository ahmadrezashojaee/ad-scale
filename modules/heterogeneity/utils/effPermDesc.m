function varargout=effPermDesc(BS,CS,TauS2, QS, Ve)
% Calculate effective permeability descriptors from streamline descriptors
% effective bulk volume and rates
% 
%
% SYNOPSIS:
% [BS, CS, TauS2]=streamtubePermDesc(GradP, LS, fluid, shortestpath);
%    
%
% PARAMETERS:
%
%   GradP        - Pressure gradients for each substep in Pollock
%                approximation
%
%   LS           - Streamline length for each substep in Pollock approximation
%
%   fluid        - MRST fluid structure where fluid.properties(1) is viscosity
%
%   shortestpath - Length of streamline if no tortuosity
%
%
% RETURNS:
%
%  BS         - Streamline/streamtube hydraulic conductance
%
%  CS         - Streamline constriction factor
%
%  TauS2      - Streamline tortuosity factor
%
%
% EXAMPLE:
%
% [BS, CS, TauS2]=streamtubePermDesc(GradP, LS, fluid, shortestpath);
%
% SEE ALSO: effPermDesc(), pollockMod()
% Written by Asgeir Nyvoll, MSc student NTNU, 2018

Q=sum(QS);
BSQS=BS.*QS;
sBSQS=sum(BSQS);
  varargout{1}=sBSQS./Ve;

if nargout>1
  varargout{2}=sum(CS.*QS)./Q;   
end

if nargout>2
  varargout{3}=sum(TauS2.*BSQS)/sBSQS;
end

if nargout>3
  varargout{4}=sum(BSQS.*(TauS2-varargout{3}(:)).^2)/sBSQS;
end

if nargout>4
  varargout{5}=sum(((1./CS-1./varargout{2}).^2).*QS)./Q;  
  

end
