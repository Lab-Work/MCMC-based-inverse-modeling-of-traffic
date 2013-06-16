%simulates the vehicle position using algorithm 2 of the article
function [updatedPostmile,vehicleSpeed] = updateVehicleWithRiemannSendingReceiving(...
    cellSpeed,...
    cellRhoCrit,...
    cellRhoMax,...
    cellWf,...
    downstreamSpeed,...
    downstreamRhoCrit,...
    downstreamRhoMax,...
    downstreamWf,...
    vmax,...
    cellIndex,...
    currentPostmile,...
    deltaX,...
    deltaT)

tol=1e-12;
debug=false;
%this is the returned variable
updatedPostmile=0;

cellRho=...
    vinverse(cellSpeed,cellRhoMax,vmax,cellWf,cellRhoCrit);
%numLanes included in rho:s

%solve flux across the cell boundary
cellSendingValue=...
    getSendingValue(...
    cellSpeed,...
    vmax,...
    cellRhoCrit,...
    cellRhoMax,...
    cellWf,1);

downstreamReceivingValue=...
    getReceivingValue(...
    downstreamSpeed,...
    vmax,...
    downstreamRhoCrit,...
    downstreamRhoMax,...
    downstreamWf,1);

cellFlux=min(downstreamReceivingValue,...
    cellSendingValue);  

%solve internal states for upstream cell
cellVcrit=vforward(cellRhoCrit,vmax,cellRhoMax,cellWf,cellRhoCrit);

cellInternalSpeed=...
    cellFlux/(cellRhoMax-cellFlux/cellWf);

cellInternalRho=vinverse(cellInternalSpeed,...
    cellRhoMax,vmax,cellWf,cellRhoCrit);


cellInternalFlow=cellInternalRho*cellInternalSpeed;

downstreamVcrit=vforward(downstreamRhoCrit,vmax,downstreamRhoMax,downstreamWf,downstreamRhoCrit);

%solve internal states for downstream cell
downstreamInternalSpeed=...
    max((-1-sqrt(1-4/vmax*cellFlux/downstreamRhoMax))/(-2/vmax),...
    (-1+sqrt(1-4/vmax*cellFlux/downstreamRhoMax))/(-2/vmax));

downstreamInternalRho=vinverse(downstreamInternalSpeed,...
    downstreamRhoMax,vmax,downstreamWf,downstreamRhoCrit);

downstreamInternalFlow=downstreamInternalRho*downstreamInternalSpeed;

downstreamRho=...
    vinverse(downstreamSpeed,downstreamRhoMax,vmax,...
    downstreamWf,downstreamRhoCrit);

downstreamFlow=downstreamRho*downstreamSpeed;
cellFlow=cellRho*cellSpeed;

%now compute the wave speeds

%interactions common to all cases
%first interaction
s1=0;                   
if(abs(cellSpeed-cellInternalSpeed)>tol)
    s1=(cellFlow-cellInternalFlow)/(cellRho-cellInternalRho);
else
    s1=0;
end
if(isnan(s1))
    keyboard
end

s2=0;
if(abs(downstreamInternalSpeed-downstreamSpeed)>tol)
    s2=(downstreamInternalFlow-downstreamFlow)/...
        (downstreamInternalRho-downstreamRho);
else
    s2=0;
end

if(isnan(s2))
    keyboard
end
s3=0;
if(abs(downstreamInternalSpeed-downstreamSpeed)>tol)
    s3=2*downstreamInternalSpeed-vmax;
else
    s3=0;
end

s4=0;
if(abs(downstreamInternalSpeed-downstreamSpeed)>tol)
    s4=2*downstreamSpeed-vmax;
else
    s4=0;
end


ta=((cellIndex+1)*deltaX-currentPostmile)/(cellSpeed-s1);
xa=(cellIndex+1)*deltaX+s1*ta;
%now just start checking into which case we fall
if(ta >= deltaT)
        updatedPostmile=currentPostmile+cellSpeed*deltaT;
        vehicleSpeed=(updatedPostmile-currentPostmile)/deltaT;
else
    tb=((cellIndex+1)*deltaX-xa+cellInternalSpeed*ta)/(cellInternalSpeed);
    if(tb<0)
        keyboard
    end
    xb=xa+cellInternalSpeed*(tb-ta);
    if(tb>=deltaT)
        updatedPostmile=xa+cellInternalSpeed*(deltaT-ta);
        vehicleSpeed=(updatedPostmile-currentPostmile)/deltaT;
        %this is final speed of the vehicle if needed
        %vehicleSpeed=cellInternalSpeed;%
        if(debug)
            keyboard
        end
    else
        if(downstreamInternalSpeed > downstreamSpeed)
            %case 1:
            tc=(downstreamInternalSpeed*tb)/...
             (downstreamInternalSpeed-s2);  
            xc=(cellIndex+1)*deltaX+s2*tc;
            if(tc>=deltaT)
                updatedPostmile=xb+downstreamInternalSpeed*(deltaT-tb);
                %this is final speed of the vehicle if needed
                %vehicleSpeed=downstreamInternalSpeed;
                vehicleSpeed=(updatedPostmile-currentPostmile)/deltaT;    
                if(debug)
                keyboard
                end
            else
                updatedPostmile=xc+downstreamSpeed*(deltaT-tc);
                %this is final speed of the vehicle if needed
                %vehicleSpeed=downstreamSpeed;
                vehicleSpeed=(updatedPostmile-currentPostmile)/deltaT;
                if(debug)
                keyboard
                end
            end
        else
            %case 2:
            tc=(downstreamInternalSpeed*tb)/...
                (downstreamInternalSpeed-s3);
            xc=(cellIndex+1)*deltaX+s3*tc;
            if(tc>=deltaT)
                updatedPostmile=xb+downstreamInternalSpeed*(deltaT-tb);
                %this is final speed of the vehicle if needed
                %vehicleSpeed=downstreamInternalSpeed;
                vehicleSpeed=(updatedPostmile-currentPostmile)/deltaT;
                if(debug)
                keyboard
                end
            else
                integrationConstant=(xc-tc*vmax-(cellIndex+1)*deltaX)...
                    /(sqrt(tc));

                td=integrationConstant^2/(s4-vmax)^2;
                xd=(cellIndex+1)*deltaX+s4*td;
                if(td>=deltaT)
                    updatedPostmile=...
                        integrationConstant*sqrt(deltaT)+...
                        deltaT*vmax+(cellIndex+1)*deltaX;
                    %this is final speed of the vehicle if needed
                    %vehicleSpeed=0.5*integrationConstant/sqrt(deltaT)+...
                    %    vmax;
                    vehicleSpeed=(updatedPostmile-currentPostmile)/deltaT;
                    if(debug)
                        keyboard
                    end
                else
                    updatedPostmile=...
                    integrationConstant*sqrt(td)+...
                        td*vmax+(cellIndex+1)*deltaX+...
                        downstreamSpeed*(deltaT-td);
                    %this is final speed of the vehicle if needed
                    %vehicleSpeed=downstreamSpeed;
                    vehicleSpeed=(updatedPostmile-currentPostmile)/deltaT;
                    if(debug)
                        keyboard
                    end
                end
            end
        end
    end
end







