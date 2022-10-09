function varargout=streamtubePermDesc(gradP, dS, fluid, s)
% Calculate streamline permeability descriptors from streamline gradients,
% line lengths, fluid model and shortest path.
% 
%
% SYNOPSIS:
% [Be,Ce,Te2,varTs2,varInvCs] = effPermDesc(BS,CS,TauS2, Qs, Ve);
%    
%
% PARAMETERS:
%
%   BS        - Streamline/streamtube hydraulic conductance
%
%   CS        - Streamline constriction factor
%
%   TauS2     - Streamline tortuosity factor
%
%   Qs        - Streamtube rate
%
%   Ve        - Effective bulk volume
%

%
% RETURNS:
%
%  Be        - Effective hydraulic conductance
%
%  Ce        - Effective constriction factor
%
%  Te2       - Effective tortuosity factor
%
%  varTs2    - Weighted variance of tortuosity factors
%
%  varInvCs  - Weighted variance of inverse constriction factors

%
% EXAMPLE:
%
% [Be,Ce,Te2,varTs2,varInvCs] = effPermDesc(BS,CS,TauS2, Qs, Ve);
%
% SEE ALSO: streamtubePermDesc()
% Written by Asgeir Nyvoll, MSc student NTNU, 2018

if length(gradP)~=length(dS)
    error('The length of pressure gradients and streamlines have to be equal')
end


gradInt=sum(1./gradP(gradP~=0).*dS(gradP~=0));

varargout{1}=-gradInt.*fluid.properties(1);

if nargout>1
    Ls=sum(dS);
    pDrop=sum(gradP.*dS);
    varargout{2}=pDrop/(Ls^2)*gradInt;
end

if nargout>2
    varargout{3}=(s/Ls)^2;
end
end