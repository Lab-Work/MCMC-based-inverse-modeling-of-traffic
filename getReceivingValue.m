function [receivingValue]=getReceivingValue(...
    velocity,vmax,rhocrit,rhomax,wf,numLanes)
%computes the receiving value of the cell
vcrit=vforward(rhocrit,vmax,rhomax,wf,rhocrit);

if (vcrit <= velocity)
    outgoingDensity=rhomax*(1-vcrit/vmax);
    receivingValue=vcrit*outgoingDensity*numLanes;
else
    outgoingDensity=rhomax*(1 / (1+velocity/wf));
    receivingValue=velocity*outgoingDensity*numLanes;
end
