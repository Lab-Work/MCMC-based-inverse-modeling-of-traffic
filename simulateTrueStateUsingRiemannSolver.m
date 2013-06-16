%this code simulates the true state 
%used in the article simulations. It models the vehicle trajectories
%using the algorithm 2 in the article
close all
clear all
%set true state
deltaT=2/3600;
vmax=77;
rhomax=[180*5,170*4,160*5];
wf=16;

vmaxForMesh=80;
deltaX=deltaT*vmaxForMesh;
domainLengthInmiles=2;
numCells=domainLengthInmiles/deltaX;

timeSteps=500
vInitial=77*ones(numCells,1);
%% plot boundary conditions
figure
vUpstream=62*ones(timeSteps+1,1);
vUpstream(50:end-300)=45;
vDownstream=15*ones(timeSteps+1,1);
vDownstream(1:40)=58;
vDownstream(220:400)=58;
timeDisc=(0:deltaT:timeSteps*deltaT);
plot(timeDisc*60,vUpstream,'k -')
hold on
 plot(timeDisc*60,vDownstream,'k --')
ylabel('\it v')
xlabel('\it t')
 %legend('upstream','downstream')
set(gca,'Ylim',[0 70])
set(gcf, 'PaperUnits', 'inches');
papersize=[3 3];
set(gcf, 'PaperSize',papersize);
width=2.5;
height=2.5;
left=(papersize(1)-width)/2;
bottom=(papersize(2)-height)/2;
myfiguresize = [left,bottom,width,height];
set(gcf, 'PaperPosition', myfiguresize);
print('-dpsc2','figs/trueBoundaryConditions.eps')
system('epstopdf figs/trueBoundaryConditions.eps')
%axis tight
%% simulate v-field with CFL=0.5
numLanes=1*ones(size(vInitial));
numLanes(15:25)=1;
numLanes=[numLanes(1);numLanes;numLanes(end)];
dropLocation=(15:25);

rhoMaxVec=zeros(numCells,1);
    rhoMaxVec(dropLocation)=rhomax(2);
    rhoMaxVec(1:dropLocation(1)-1)=rhomax(1);
   rhoMaxVec(dropLocation(end)+1:end)=rhomax(3);
rhoMaxVec=[rhoMaxVec(1);rhoMaxVec;rhoMaxVec(end)];

rhoCritVec=rhoMaxVec.*(wf/vmax);

vupdated=updatevHalfCFL(vInitial,deltaX,deltaT,timeSteps,vDownstream,vUpstream,rhoCritVec,vmax,wf*ones(numCells+2,1),rhoMaxVec,numLanes);
%% plotting
cmap=flipud(colormap(jet));
spaceDisc=0:deltaX:(domainLengthInmiles)-deltaX;
timeDisc=(0:deltaT:timeSteps*deltaT);
figure
imagesc(timeDisc*60,spaceDisc,vupdated)
colormap(cmap)
set(gca,'Clim',[0 80])
colorbar('location','south')
axis xy
axis tight
xlabel('{\it t}')
ylabel('postmile')
set(gcf, 'PaperUnits', 'inches');
papersize=[3 3];
set(gcf, 'PaperSize',papersize);
width=2.5;
height=2.5;
left=(papersize(1)-width)/2;
bottom=(papersize(2)-height)/2;
myfiguresize = [left,bottom,width,height];
set(gcf, 'PaperPosition', myfiguresize);
print('-dpsc2','figs/trueState.eps')
system('epstopdf figs/trueState.eps')

%% now we simulate the vehicle trajectories
%we do this one vehicle at a time to keep the code simple
% finally we plot the trajectories
figure
X=[];
V=[];
cellNumbers=[];
for ii=1:14:length(timeDisc)
    entryTime=(ii-1)*deltaT+0.00000001;
    entryPostmile=0.000001;
    [vupdated,v,x,t,cellNumber]=updatevHalfCFLAndSimulateVehicle(vInitial...
        ,deltaX...
        ,deltaT...
        ,timeSteps...
        ,vDownstream...
        ,vUpstream...
        ,rhoCritVec...
        ,vmax...
        ,wf*ones(numCells+2,1)...
        ,rhoMaxVec...
        ,numLanes...
        ,entryTime);
    X=[X,x];
    V=[V,v];
    cellNumbers=[cellNumbers,cellNumber];
    hold on
    plot(t*60,x,'-','color',[0.4 0.4 0.4])
end
axis tight
%grid on
set(gcf, 'PaperUnits', 'inches');
papersize=[3 3];
set(gcf, 'PaperSize',papersize);
width=2.5;
height=2.5;
left=(papersize(1)-width)/2;
bottom=(papersize(2)-height)/2;
myfiguresize = [left,bottom,width,height];
set(gcf, 'PaperPosition', myfiguresize);
xlabel('{\it t}')
ylabel('postmile')

print('-dpsc2','figs/trueTrajectories.eps')
system('epstopdf figs/trueTrajectories.eps')
%%
vtrue=vupdated;
save vtrueRiemann vtrue timeDisc spaceDisc
Vprobes=V;
Xprobes=X;
XprobesRiemann=X;
VprobesRiemann=V;
save probemeasurementsRiemann Xprobes Vprobes t vtrue cellNumbers


