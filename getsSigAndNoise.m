function [theSig allCellsNoise theAct numBins numReps] = getsSigAndNoise(thisNeu);
%%%% thisNeu is the spike count for all trials (columns=bins, rows=events)
numBins=size(thisNeu,1);  %%%% total number of bins
numReps=size(thisNeu,2); %%%% number of repetitions of the movie
theSig=mean(thisNeu')';  %%% extract the signal by taking the mean
noiseAux=demeaner(thisNeu'); %%% extracts the signal from the total activity to calculate the noise
allCellsNoise=noiseAux(:);  %%% concatenates the noise
theAct=thisNeu(:);  %%% concatenates the activity