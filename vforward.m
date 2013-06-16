function [v] = vforward(rho,vmax,rhomax,wf,rhocrit)
%VFORWARD Evaluate the mapping from rho space to v space
%rhocrit=vinverse(vcrit,vmax,rhomax,wf,vcrit);
v=0;
if rho<= rhocrit
    v=vmax*(1-rho/rhomax);
elseif rho>rhocrit && rho<=rhomax
    v=-wf*(1-rhomax/rho);
elseif rho<0
    disp('Warning rho < 0 in vforward');
    v=vmax;
end
    