function [probeDraws]=getProbeDraws(lowerBound,probeSigma,meanValue,numSamples)

probeDraws=randn(1,numSamples)*probeSigma+meanValue;
samplesNotGood=true;
while(samplesNotGood)
        if(~isempty(find(probeDraws<lowerBound)))
            samplesNotGood=true;
            probeDraws=meanValue+probeSigma*...
                randn(1,numSamples);
        else
            samplesNotGood=false;
            
        end
end