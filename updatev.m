%simulates velocity field using godunov scheme
%all parameter vector inputs (rhocrit, rhomax, vmax, wf, numlanes) assume ghost cells in the first and last cell
%except the initial state v
function vupdated=updatev(v,deltaX,deltaT,timeSteps,vDownstream,vUpstream,rhocrit,vmax,wf,rhomax,numLanes)

numberOfActualCells=length(v);
vupdated=zeros(length(v)+2,timeSteps+1);
vupdated(2:end-1,1)=v;
vupdated(1,:)=vUpstream';
vupdated(end,:)=vDownstream';

for tt=2:timeSteps+1
    
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

vupdated=vupdated(2:end-1,:);
