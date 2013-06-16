
function vupdated=updatevHalfCFL(v,deltaX,deltaT,timeSteps,vDownstream,vUpstream,rhocrit,vmax,wf,rhomax,numLanes)
%%updates v field with half of the cfl specified
%by deltat and delta x
%uses linear interpolation for boundary conditions
numberOfActualCells=length(v);
vupdated=zeros(length(v)+2,2*(timeSteps)+1);
vupdated(2:end-1,1)=v;

vupdated(1,:)=interp1(1:(timeSteps+1),vUpstream',1:0.5:(timeSteps+1));
vupdated(end,:)=interp1(1:(timeSteps+1),vDownstream',1:0.5:(timeSteps+1));
%start godunov scheme
for tt=2:2*(timeSteps)+1
    
    for cell=2:numberOfActualCells+1
        
        upstreamFlux = getSendingValue(vupdated(cell-1,tt-1),vmax,rhocrit(cell-1),rhomax(cell-1),wf(cell-1),numLanes(cell-1));
        upstreamFlux=min(upstreamFlux,...
            getReceivingValue(vupdated(cell,tt-1),vmax,rhocrit(cell),rhomax(cell),wf(cell),numLanes(cell)));
        
        downstreamFlux = getSendingValue(vupdated(cell,tt-1),vmax,rhocrit(cell),rhomax(cell),wf(cell),numLanes(cell));
        
        downstreamFlux=min(downstreamFlux,getReceivingValue(vupdated(cell+1,tt-1),vmax,rhocrit(cell+1),rhomax(cell+1),wf(cell+1),numLanes(cell+1)));
        
        updatedrho=vinverse(vupdated(cell,tt-1),rhomax(cell),vmax,wf(cell),rhocrit(cell))-...
            deltaT/deltaX*(...
            (downstreamFlux-upstreamFlux)/numLanes(cell));    

        vupdated(cell,tt)=vforward(...
             updatedrho,vmax,rhomax(cell),wf(cell),rhocrit(cell));
    
    end

end

vupdated=vupdated(2:end-1,1:2:end);