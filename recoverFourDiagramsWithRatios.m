%THIS IS THE RUN THAT RECOVERS THE FOUR DIAGRAMS
%IN THE ARTICLE: three ratios + one jamdensity + wf + vmax
clear all
close all
noiseStd=2
load vtrueRiemann
load probemeasurementsRiemann
Vprobes=Vprobes+noiseStd*randn(size(Vprobes));
save VprobesNoisy Vprobes
load VprobesNoisy

%set parameters
deltaT=2/3600;
timeSteps=500
vUpstream=62*ones(timeSteps+1,1);
vUpstream(50:end-300)=45;
vDownstream=15*ones(timeSteps+1,1);
vDownstream(1:40)=58;
vDownstream(220:400)=58;
vmaxMesh=80;
base=150*5;
rhomax=[base;base;base;base+75]
wf=18;
vmax=75;
deltaX=deltaT*vmaxMesh;%1/16;
domainLengthInmiles=2;
numCells=domainLengthInmiles/deltaX;
rhocrit=rhomax.*(wf/vmax);

%initial state

vInitial=77*ones(numCells,1);

numCells=domainLengthInmiles/deltaX;

burnIn=100000;
samples=200000

upperRho=250;
lowerRho=120;
%mixing coefficients are set here
mixingcoefficients=[0.1,0.05];
Rhocrit=zeros(samples,4);
Vinit=zeros(samples,1);
Wf=zeros(samples,1);

numLanes=ones(numCells+2,1);
%set rhomax values into right cells
dropLocation=(15:25);
rhoMaxVec=zeros(numCells,1);
rhoMaxVec(dropLocation)=rhomax(2);
rhoMaxVec(1:dropLocation(1)-1)=rhomax(1);
rhoMaxVec(dropLocation(end)+1:end-10)=rhomax(3);
rhoMaxVec(end-9:end)=rhomax(4);
rhoMaxVec=[rhoMaxVec(1);rhoMaxVec;rhoMaxVec(end)];
rhoCritVec=rhoMaxVec.*(wf/vmax);

wfVec=wf*ones(numCells+2,1);

%get first v-field with initial condition
vupdated=updatev(vInitial,deltaX,deltaT,timeSteps,vDownstream,vUpstream,rhoCritVec,vmax,wfVec,rhoMaxVec,numLanes);
%figure,imagesc(vupdated-vtrue),colorbar
%axis xy
newMeasurements=[];
probeMeasurements=[];
for timeStep=1:timeSteps+1
    %find the cells that have vehicles at timeStep
    vehicleIndices=~isnan(cellNumbers(timeStep,:));
    cellsWithVehicles=cellNumbers(timeStep,vehicleIndices);
    newMeasurements=[newMeasurements;vupdated(cellsWithVehicles,timeStep)];
    probeMeasurements=[probeMeasurements;Vprobes(timeStep,vehicleIndices)'];
end

noiseStd=3.0;
noiseDiagonals=noiseStd^2*ones(size(newMeasurements));
%compute log likelihood + posterior
logLikeOld=-0.5*(newMeasurements-probeMeasurements)'*((newMeasurements-probeMeasurements)./noiseDiagonals);
logLikeOld=logLikeOld+log(evaluateUniformDistribution(rhomax(2),110*4,250*5));
likelihoods=zeros(samples,1);
acceptedSamples=0;
ratio=[rhomax(2)/rhomax(1);rhomax(3)/rhomax(1);rhomax(4)/rhomax(3)];
%main loop
for sample=1:samples
    accepted=false;
    
    rhomaxproposal=zeros(4,1);
    rhomaxproposal(2)=rhomax(2)+randn(1,1)*10;
    ratioproposal=ratio+randn(3,1)*0.000075;
    while(~isempty(find(ratioproposal <0)))
        ratioproposal=ratio+randn(2,1)*0.000075;
    end
    rhomaxproposal(1)=rhomaxproposal(2)*ratioproposal(1);
    rhomaxproposal(3)=rhomaxproposal(2)*ratioproposal(2);
    rhomaxproposal(4)=rhomaxproposal(3)*ratioproposal(3);
    
    wfproposal=wf+randn(1)*mixingcoefficients(2);
    while(wfproposal <5 || wfproposal > 25)
        wfproposal=wf+randn(1)*mixingcoefficients(2);
    end
    
    vmaxproposal=vmax+randn(1)*mixingcoefficients(1);
    while(vmaxproposal <0.1 || vmaxproposal > vmaxMesh)
        vmaxproposal=vmax+randn(1)*mixingcoefficients(1);
    end
    rhoMaxVec=zeros(numCells,1);
    rhoMaxVec(dropLocation)=rhomaxproposal(2);
    rhoMaxVec(1:dropLocation(1)-1)=rhomaxproposal(1);
   rhoMaxVec(dropLocation(end)+1:end-10)=rhomaxproposal(3);
    rhoMaxVec(end-9:end)=rhomaxproposal(4);
    rhoMaxVec=[rhoMaxVec(1);rhoMaxVec;rhoMaxVec(end)];
    rhoCritVec=rhoMaxVec.*(wfproposal/vmaxproposal);
    
    wfVec=wfproposal*ones(numCells+2,1);
    vupdated=updatev(vInitial,deltaX,deltaT,timeSteps,vDownstream,vUpstream,rhoCritVec,vmaxproposal,wfVec,rhoMaxVec,numLanes);
    
    newMeasurements=[];
    probeMeasurements=[];
    for timeStep=1:timeSteps+1
        %find the cells that have vehicles at timeStep
        vehicleIndices=~isnan(cellNumbers(timeStep,:));
        cellsWithVehicles=cellNumbers(timeStep,vehicleIndices);
        newMeasurements=[newMeasurements;vupdated(cellsWithVehicles,timeStep)];
        probeMeasurements=[probeMeasurements;Vprobes(timeStep,vehicleIndices)'];
    end
    logLikeNew=-0.5*(newMeasurements-probeMeasurements)'*((newMeasurements-probeMeasurements)./noiseDiagonals);
    
    logPrior=log(evaluateUniformDistribution(rhomaxproposal(2),110*4,250*5));
    logLikeNew=logLikeNew+logPrior;
    alpha=exp(logLikeNew-logLikeOld);
    %check if accepted
    if (rand(1)<alpha)
        accepted=true;
        vmax=vmaxproposal;
        rhomax=rhomaxproposal;
        ratio=ratioproposal;
        wf=wfproposal;
        logLikeOld=logLikeNew;
    end
    
    Rhocrit(sample,:)=rhomax;
    Wf(sample)=wf;
    Vinit(sample)=vmax;
    likelihoods(sample)=logLikeOld;
    if(accepted==true)
        acceptedSamples=acceptedSamples+1;
    end
    
    if (sample>0 && mod(sample,100)==0)
        disp(['Round: ' num2str(sample)])
        disp(['Acceptance ratio: ' num2str(acceptedSamples/sample*100)])
        
    end
    
end
save results4Diagrams Rhocrit Vinit Wf likelihoods burnIn
%%

figure
cmap=flipud(colormap(jet));
spaceDisc=0:deltaX:(domainLengthInmiles)-deltaX;
timeDisc=(0:deltaT:timeSteps*deltaT);
imagesc(timeDisc*60,spaceDisc,vupdated)
axis xy
axis tight
colormap(cmap)
set(gca,'Clim',[0 80])
%% plotting without running the simulation again works also
close
load results4Diagrams
figure
subplot(2,2,1)
hold on 
box on
numSamples=length(Rhocrit(:,2));


truerho=[180*5,170*4,160*5];
plot([1 numSamples],truerho(2)*ones(1,2),'-','Color',[0.7 0.7 0.7],'linewidth',3)
plot(Rhocrit(:,1),'k -')

ylabel('\rho_{max}')
xlabel('sample')
subplot(2,2,[2])
hold on
box on
numSamples=length(Rhocrit(:,1));

truerho=[180*5,170*4,160*5,160*5];
plot([1 numSamples],truerho(1)/truerho(2)*ones(1,2),'-','Color',[0.7 0.7 0.7],'linewidth',3)
plot([1 numSamples],truerho(3)/truerho(2)*ones(1,2),'--','color',[0.7 0.7 0.7],'linewidth',3)
plot([1 numSamples],truerho(4)/truerho(3)*ones(1,2),'-.','color',[0.7 0.7 0.7],'linewidth',3)

trueratio1=truerho(1)/truerho(2)
trueratio2=truerho(3)/truerho(2)
trueratio3=truerho(4)/truerho(3)
plot(Rhocrit(:,1)./Rhocrit(:,2),'k -')
plot(Rhocrit(:,3)./Rhocrit(:,2),'k --')
plot(Rhocrit(:,4)./Rhocrit(:,3),'k -.')
ratio1=Rhocrit(:,1)./Rhocrit(:,2);
ratio1mean=mean(ratio1(burnIn:end))
ratio1std=std(ratio1(burnIn:end))

ratio2=Rhocrit(:,3)./Rhocrit(:,2);
ratio2mean=mean(ratio2(burnIn:end))
ratio2std=std(ratio2(burnIn:end))
ratio3=Rhocrit(:,4)./Rhocrit(:,3);
ratio3mean=mean(ratio3(burnIn:end))
ratio3std=std(ratio3(burnIn:end))

rhomax2mean=mean(Rhocrit(burnIn:end,2))
rhomax2std=std(Rhocrit(burnIn:end,2))

ylabel('\rho_{max} ratio')
xlabel('sample')
subplot(2,2,3)
wftrue=16;
hold on
box on
plot([1 numSamples],wftrue*[1 1],'-','Color',[0.7 0.7 0.7],'linewidth',3)

plot(Wf,'k -')
ylabel('\it w_f')
xlabel('sample')
subplot(2,2,4)
hold on
box on
vmaxtrue=77;
plot([1 numSamples],vmaxtrue*[1 1],'-','Color',[0.7 0.7 0.7],'linewidth',3)

plot(Vinit,'k -')
ylabel('{\it v}_{max}')
xlabel('sample')
rhocritmean=mean(Rhocrit(burnIn:end,:))
rhocritstd=std(Rhocrit(burnIn:end,:))
wfmean=mean(Wf(burnIn:end))
wfstd=std(Wf(burnIn:end))
vmaxmean=mean(Vinit(burnIn:end))
vmaxstd=std(Vinit(burnIn:end))
set(gcf, 'PaperUnits', 'inches');
papersize=[5 5];
set(gcf, 'PaperSize',papersize);
width=4.5;
height=4.5;
left=(papersize(1)-width)/2;
bottom=(papersize(2)-height)/2;
myfiguresize = [left,bottom,width,height];
set(gcf, 'PaperPosition', myfiguresize);
print('-dpsc2','fourDiagramChains.eps')
system('epstopdf fourDiagramChains.eps')


