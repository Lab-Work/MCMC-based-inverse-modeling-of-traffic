function [laneVector]=getLaneVector(numberOfMainLanes,numberOfAnomalyLanes,startPoint,endPoint,numberOfCells,domainLength)

laneVector=numberOfMainLanes*ones(numberOfCells,1);

startIndex=floor(startPoint/domainLength*numberOfCells);

endIndex=round(endPoint/domainLength*numberOfCells);


laneVector(startIndex:endIndex)=numberOfAnomalyLanes;

