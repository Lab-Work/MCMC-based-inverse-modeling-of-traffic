function [rho] = vinverse(v,rhomax,vmax,wf,rhocrit)
%VFORWARD Evaluate the mapping from v space to rho space
vcrit=vforward(rhocrit,vmax,rhomax,wf,rhocrit);
rho=0;

if v>= vcrit && v<vmax
    rho=rhomax*(1-v/vmax);
elseif v < vcrit && v >=0
    rho=rhomax*(1/( 1+v/wf ));
elseif v<0
    disp('Warning v<0 in vinverse');
    rho=rhomax;
end