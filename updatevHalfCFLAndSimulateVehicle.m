
function [vupdated,v,x,t,cellNumbers]=updatevHalfCFLAndSimulateVehicle(v,deltaX,deltaT,timeSteps,vDownstream,vUpstream,rhocrit,vmax,wf,rhomax,numLanes,entryTime)
%updates v field with half of the cfl specified by deltat and delta x
%uses linear interpolation for boundary conditions
%also simulates the vehicle trajectory from a given initial condition (x,t)
%returns the full velocity field with speed vector of the vehicle for
%each time step, as well as position for each time step, time
%discretization used, and cell number of the vehicle for each time step.
%NaN means vehicle was not in the domain at that time step
numberOfActualCells=length(v);

%these lines are for cfl=0.5
vupdated=zeros(length(v)+2,2*(timeSteps)+1);
deltaT=deltaT/2;
%timeSteps=timeSteps*2;
maxTimeIndex=timeSteps;

%this line is for cfl=1
%vupdated=zeros(length(v)+2,timeSteps+1);
%maxTimeIndex=timeSteps;

vupdated(2:end-1,1)=v;
%set parameters for vehicle simulation 
entryPostmile=0.000001;
spaceDisc=0:deltaX:(deltaX*numberOfActualCells);
currentCellNumber=intersect(find(spaceDisc<=entryPostmile),...
        find(spaceDisc>=entryPostmile-deltaX));
currentPostmile=entryPostmile;%spaceDisc(currentCellNumber);

maxCellNumber=size(v,1);


x=zeros(maxTimeIndex*2+1,1)*NaN;
v=zeros(maxTimeIndex*2+1,1)*NaN;
cellNumbers=zeros(maxTimeIndex*2+1,1)*NaN;
t=0:deltaT:((maxTimeIndex*2)*deltaT);
timeOfEntryIndex=intersect(find(t<=entryTime),find(t>=entryTime-deltaT));
if(isempty(currentCellNumber) || isempty(timeOfEntryIndex))
    keyboard;
end

%these lines for cfl =0.5
vupdated(1,:)=interp1(1:(timeSteps+1),vUpstream',1:0.5:(timeSteps+1));
vupdated(end,:)=interp1(1:(timeSteps+1),vDownstream',1:0.5:(timeSteps+1));

%these lines for cfl=1
%vupdated(1,:)=vUpstream';
%vupdated(end,:)=vDownstream';

for tt=2:2*timeSteps+1
    
    updatedForTimeStep=false;
    for cell=2:numberOfActualCells+1
        
        updatePostmile=false;
        %current cell number refes to index without ghost cells
        if (currentCellNumber==cell-1 && tt>=timeOfEntryIndex && ~updatedForTimeStep)
            updatePostmile=true;
        end
        if (currentPostmile>spaceDisc(end))
            updatePostmile=false;
        end
        
        if (updatePostmile)
            %here goes the postmile update
            %using riemann problem formulation
            
            [updatedPostmile,vehicleSpeed]=updateVehicleWithRiemannSendingReceiving(...
                vupdated(cell,tt-1),...
                rhocrit(cell),...
                rhomax(cell),...
                wf(cell),...                       
                vupdated(cell+1,tt-1),...
                rhocrit(cell+1),...
                rhomax(cell+1),...
                wf(cell+1),...
                vmax,...
                cell-2,...
                currentPostmile,...
                deltaX,...
                deltaT);
            currentPostmile=updatedPostmile;
            if (abs(currentPostmile)>real(currentPostmile))
                keyboard
            end
            
            if (currentPostmile < spaceDisc(end))
            
                x(tt-1)=currentPostmile;
            %find the current cell number again
    
             currentCellNumber=intersect(find(spaceDisc<=currentPostmile),...
                find(spaceDisc>=currentPostmile-deltaX));
                if (length(currentCellNumber) ==0)
                    keyboard
                end

                cellNumbers(tt-1)=currentCellNumber;
           
             v(tt-1)=vehicleSpeed;
            end
            updatedForTimeStep=true;
            %keyboard
        end
        
        upstreamFlux = getSendingValue(vupdated(cell-1,tt-1),vmax,rhocrit(cell-1),rhomax(cell-1),wf(cell-1),numLanes(cell-1));
        
        upstreamFlux=min(upstreamFlux,...
            getReceivingValue(vupdated(cell,tt-1),vmax,rhocrit(cell),rhomax(cell),wf(cell),numLanes(cell)));
        
        downstreamFlux = getSendingValue(vupdated(cell,tt-1),vmax,rhocrit(cell),rhomax(cell),wf(cell),numLanes(cell));
        
        downstreamFlux=min(downstreamFlux,...
            getReceivingValue(vupdated(cell+1,tt-1),vmax,rhocrit(cell+1),rhomax(cell+1),wf(cell+1),numLanes(cell+1)));
        
        updatedrho=vinverse(vupdated(cell,tt-1),rhomax(cell),vmax,wf(cell),rhocrit(cell))-...
            deltaT/deltaX*(...
            (downstreamFlux-upstreamFlux)/numLanes(cell));    

        vupdated(cell,tt)=vforward(...
             updatedrho,vmax,rhomax(cell),wf(cell),rhocrit(cell));
    
    end

end

%this line is for cfl=1;
%vupdated=vupdated(2:end-1,:);
%this line is for cfl=0.5
vupdated=vupdated(2:end-1,[1,3:2:end]);

v=v([1,3:2:end]);
x=x([1,3:2:end]);

t=t([1,3:2:end]);

cellNumbers=cellNumbers([1,3:2:end]);
