% clear all
% load windowExpsTreatCont %%%% loads perievent histograms
% windowTreatOrCont=windowExpsTreat;
clear all
load theExperiments %%%% loads perievent histograms in variable experiments
%%%% theExperiemnts contains spike counts for each different bin size. It
%%%% contains information for each experiment, for each neuron, for each
%%%% type of stimulus, for each block repetition (3 blocks that are further
%%%% concatenated), for each stimulus repetition (30 of them), and for each
%%%% experimental condition (nb stimulation vs control). Look below for the
%%%% structure of the variable.
nmbrBtstr=10;%% number of bootstraps (keep low for quick runs of the code, increase for accuracy)
for binSizeNum=1:length(theExperiments)-1 %%% goes through each window size
    binSize(binSizeNum)=theExperiments{binSizeNum+1}.window;  %%% windowSize
    experiments=theExperiments{binSizeNum+1}.exp;  %%%% binning at that windowSize  
    %%%%%%%%%%%%%%%%%%%%%%%% gets covariance matrices %%%%%%%%%%%%%%%%
    for expNum=1:length(experiments) %%%% loops through all the experiments
        for nbStim=1:2; %%% 1 is Ach, 2 is control
            clear varSig,clear allCellsNoise;clear theSig;clear theAct; clear allSigsAux
            for cellNum=1:size(experiments{expNum},2) %%% goes through each neuron
                thisCell=[experiments{expNum}{cellNum}.periEvent{nbStim,1}...
                    experiments{expNum}{cellNum}.periEvent{nbStim,2}...
                    experiments{expNum}{cellNum}.periEvent{nbStim,3}]; %%%% concatenates the three blocks of 5 seconds each
                [theSigAux allCellsNoiseAux theActAux numBins numReps]=getsSigAndNoise(thisCell');%%% extracts signal and noise
                theSig(:,cellNum)=theSigAux; %%% all signals
                allCellsNoise(:,cellNum)=allCellsNoiseAux; %%% all noises (concatenated for each neuron)
                theAct(:,cellNum)=theActAux;
                meanSpikesAux(cellNum)=mean(theActAux);
            end
            dfNoise=(numBins*(numReps-1));  %%%%%% degrees of freedom of the noise
            covNoise{expNum}{nbStim}=(numBins*numReps)*cov(allCellsNoise,1)/dfNoise;%%% unbiased estimator of the noise covariance matrix
            covSig{expNum}{nbStim}=cov(theSig);
            [corrSig{expNum}{nbStim} stdSig{expNum}{nbStim}]=corrcov(covSig{expNum}{nbStim});%%% signal correlations
            [corrNoise{expNum}{nbStim} stdNoise{expNum}{nbStim}]=corrcov(covNoise{expNum}{nbStim});%%%% noise correlations
            meanSpikes{expNum}{nbStim}=meanSpikesAux;            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    allSigsCorrs{1}=[];allSigsCorrs{2}=[];allNoiseCorrs{1}=[];allNoiseCorrs{2}=[];
    minSigNoise{1}=[];minSigNoise{2}=[];expID=[];allSigs{1}=[];allSigs{2}=[];allNoises{1}=[];allNoises{2}=[];
    allRates{1}=[];allRates{2}=[];
    for expNum=[1:length(experiments)] %%% loops through experiments
        stdSigNB=stdSig{expNum};
        stdNoiseNB=stdNoise{expNum};
        corrSigNB=corrSig{expNum};
        corrNoiseNB=corrNoise{expNum};
        
        for nbStim=1:2      %%% 1 is Ach, 2 is control
            corrSigAux=offLowerTri(corrSigNB{nbStim})';
            corrNoiseAux=offLowerTri(corrNoiseNB{nbStim})';
            allSigsCorrs{nbStim}=[allSigsCorrs{nbStim} corrSigAux];
            allNoiseCorrs{nbStim}=[allNoiseCorrs{nbStim} corrNoiseAux];
            allSigs{nbStim}=[allSigs{nbStim}; stdSigNB{nbStim}];
            allNoises{nbStim}=[allNoises{nbStim}; stdNoiseNB{nbStim}];
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%% CALCULATES VARIABLES IN THE PAPER %%%%%%%%%%%%
    %%% percentile change in signal amplitude
    percSigs=allSigs{1}./allSigs{2}*100-100; 
    sigBinProp(binSizeNum)=median(percSigs);
    sigBinPropConf(binSizeNum,:)=prctile(bootstrp(nmbrBtstr,@median,percSigs),[2.5 97.5])-median(percSigs);
    %%% percentile change in noise amplitude
    percNoise=allNoises{1}./allNoises{2}*100-100; 
    noiseBinProp(binSizeNum)=median(percNoise);
    noiseBinPropBtstr=bootstrp(nmbrBtstr,@median,percNoise);
    noiseBinPropConf(binSizeNum,:)=prctile(noiseBinPropBtstr,[2.5 97.5])-median(percNoise);
    %%% change in signal correlation
    diffSigCorrsAux=allSigsCorrs{1}-allSigsCorrs{2}; 
    diffSigCorrs(binSizeNum)=median(diffSigCorrsAux);
    diffSigCorrsBtstr=bootstrp(nmbrBtstr,@median,diffSigCorrsAux);
    diffSigCorrsConf(binSizeNum,:)=prctile(diffSigCorrsBtstr,[2.5 97.5])-median(diffSigCorrsAux);
    %%% change in noise correlation
    diffNoiseCorrsAux=allNoiseCorrs{1}-allNoiseCorrs{2};  
    diffNoiseCorrs(binSizeNum)=median(diffNoiseCorrsAux);
    diffNoiseCorrsBtstr=bootstrp(nmbrBtstr,@median,diffNoiseCorrsAux);
    diffNoiseCorrsConf(binSizeNum,:)=prctile(diffNoiseCorrsBtstr,[2.5 97.5])-median(diffNoiseCorrsAux);
    %%% correlations slope 
    [slopeAch(binSizeNum), slopeCon(binSizeNum),confSlope1(binSizeNum,:),confSlope2(binSizeNum,:),slopePvalNum(binSizeNum),slopeAchStd(binSizeNum),slopeConStd(binSizeNum)]=bootStrapDemingConfPNAS(allSigsCorrs{1},allNoiseCorrs{1},allSigsCorrs{2},allNoiseCorrs{2},100)
   
end
%

%%%%%%%%%%%%%% FOR DIFFERENCES IN SIGNAL AND NOISE %%%%%%%%%%%%%%%
figure
scatter(binSize,sigBinProp,150,'k','fill'),hold on,scatter(binSize,noiseBinProp,150,'k')
errorbar(binSize,sigBinProp,sigBinPropConf(:,1),sigBinPropConf(:,2),'ko','MarkerSize',10),
errorbar(binSize,noiseBinProp,noiseBinPropConf(:,1),noiseBinPropConf(:,2),'.k','MarkerSize',10)
plot([0 1.4],[0 0],'k')
% plotsave('amplitudeIncreaseErrorBars')

%%%%%%%%%%%% FOR PLAIN CHANGE IN SIGNAL AND NOISE CORRELATIONS %%%%%%%
figure
scatter(binSize,diffSigCorrs,150,'k','fill'),hold on,scatter(binSize,diffNoiseCorrs,150,'k')
errorbar(binSize,diffSigCorrs,diffSigCorrsConf(:,1),diffSigCorrsConf(:,2),'ko','MarkerSize',10),
errorbar(binSize,diffNoiseCorrs,diffNoiseCorrsConf(:,1),diffNoiseCorrsConf(:,2),'.k','MarkerSize',10)
plot([0 1.4],[0 0],'k')

%%%%%%%%%% FOR SLOPES AND BOOTSTRAP CONFIDENCE INTERVALS %%%%%%%%%%%%%%%
figure
scatter(binSize,slopeCon,150,'k'),hold on,scatter(binSize,slopeAch,150,'k','fill')
errorbar(binSize,slopeCon,confSlope2(:,1),confSlope2(:,2),'k.'),
errorbar(binSize,slopeAch,confSlope1(:,1),confSlope1(:,2),'k.')









