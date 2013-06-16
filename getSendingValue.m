function [sendingValue]=getSendingValue(...
    velocity,vmax,rhocrit,rhomax,wf,numLanes)
%computes the sending value of the cell
 vcrit=vforward(rhocrit,vmax,rhomax,wf,rhocrit);

if (vcrit <= velocity)
   incomingDensity=rhomax*(1-velocity/vmax);
   sendingValue=incomingDensity*velocity*numLanes;
else
   
    incomingDensity=rhomax*( 1 / (1+vcrit/wf) );
    sendingValue=vcrit*incomingDensity*numLanes;
end