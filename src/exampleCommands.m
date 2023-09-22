
cd('/Users/RobertBeardmore/Dropbox/EUCAST-ECOFF')

%%

clc
clear all
close all

load('./mat/all-EUCAST-Histograms.mat')
load('./mat/extractedEUCAST-SIR-Data.mat')
load('./mat/EUCASTCLSImatches.mat');
load('./neuralNets/bestNets.mat')
load('./mat/bothBreakpointSets.mat');
load('./mat/allnewHistosANNdecisions.mat')
load('./mat/ANNdecisions.mat')
rawCLSIbreakpoints = importCLSIbreakpoint('./CLSI/breakpoints.csv');

showbreakpointBreakdown;

%%

ECOFFexamples();

%%

%this is out of date and superceded by the next section
%plotCLSIEUCASTcorrelation(consolidatedCLSIEUCASTMatches);

%%

plotCLSIEUCASTcorrelation2(bothBreakpointSets);

%%

summariseSIRdata(consolidatedSIRData);

%%

spectralEmbedding(allEUCASTHistograms,consolidatedSIRData);

%%

ignoreList = analyseSIR1bump(consolidatedSIRData);
analyseSIR1bump(consolidatedSIRData,ignoreList);

%%

analyseOneHistogram(allEUCASTHistograms,'Ampicillin','Haemophilus influenzae ATCC 49247')
analyseOneHistogram(allEUCASTHistograms,'Cefoxitin','Staphylococcus aureus MRSA')
analyseOneHistogram(allEUCASTHistograms,'Ceftazidime','Providencia stuartii')
analyseOneHistogram(allEUCASTHistograms,'Cefotaxime','Klebsiella aerogenes')
analyseOneHistogram(allEUCASTHistograms,'Ceftolozane','Citrobacter freundii')
analyseOneHistogram(allEUCASTHistograms,'Cefiderocol','Escherichia coli')

%%

analyseSIR2bump(consolidatedSIRData,1);
analyseSIR2bump(consolidatedSIRData,2);
analyseSIR2bump(consolidatedSIRData,3);
analyseSIR2bump(consolidatedSIRData,4);
analyseSIR2bump(consolidatedSIRData,5);

%%

analyseSIR3bump(consolidatedSIRData,1);
[ignoreListOut,mdl,tbl,targetString] = analyseSIR3bump(consolidatedSIRData,2);
[ignoreListOut,mdl,tbl,targetString] = analyseSIR3bump(consolidatedSIRData,3);

%%

regressionLearner(tbl,targetString)

%%

defaultNNensemble = trainNNonAllHistogramData();

%% CLSI ANNE training
ensembleSize = 100;
for J = 1:5
    yourname = ['rob',num2str(J)];
    
    poolobj = gcp;
    
    for NNtype = 1:6
        NNensemble = trainNNonAllHistogramData(NNtype,ensembleSize);
        save(['./neuralNets/c_',yourname,'_NNensembleNtype',num2str(NNtype),'.mat'],'NNensemble');
    end
end

%% EUCAST ANNE training
ensembleSize = 100;
for J = 1:5
    yourname = ['rob',num2str(J)];
    
    poolobj = gcp;
    
    for NNtype = 1:6
        NNensemble = trainNNonAllHistogramData(NNtype,ensembleSize);
        save(['./neuralNets/e_',yourname,'_NNensembleNtype',num2str(NNtype),'.mat'],'NNensemble');
    end
end

%%

testBestNNs(bestNets,lDosages,1);
testBestNNs(bestNets,lDosages,2);
testBestNNs(bestNets,lDosages,3);

%%

% Shows correlations of ANN BP S and R decisions against CLSI ones
% Linear corrected ANN decisions based on CLSI data are shown here, but
% currently are not used later to generate data in newHistosANNdecisionData
% (for histograms where CLSI have not already made decisions, likely because
% drug-bug pairs are not in clinical usage).

compareANNwithCLSIbp(ANNdecisionData,1,'r');
compareANNwithCLSIbp(newHistosANNdecisionData,3,[1 1 1]/2);

%%

%These are currently the same figure. They are only different if some form
%of validated post-ANN-decision correction is used, but this currently is OFF.

ANNanalyseOneHistogram('Cefiderocol','Escherichia coli',newHistosANNdecisionData,2)

%%

ANNanalyseOneHistogram('Vancomycin','Staphylococcus aureus',ANNdecisionData)
ANNanalyseOneHistogram('Ceftazidime','Klebsiella pneumoniae',ANNdecisionData)
ANNanalyseOneHistogram('Ampicillin','Klebsiella pneumoniae',ANNdecisionData)
ANNanalyseOneHistogram('Tigecycline','Klebsiella pneumoniae',ANNdecisionData)
ANNanalyseOneHistogram('Tigecycline','Staphylococcus aureus',ANNdecisionData)
ANNanalyseOneHistogram('Tigecycline','Escherichia coli',ANNdecisionData)
ANNanalyseOneHistogram('Ertapenem','Escherichia coli',ANNdecisionData)
ANNanalyseOneHistogram('Ertapenem','Klebsiella pneumoniae',ANNdecisionData)

%%

%where no CLSI decision exists:
ANNanalyseOneHistogram('Tetracycline','Escherichia coli',newHistosANNdecisionData)
ANNanalyseOneHistogram('Amoxicillin','Escherichia coli',newHistosANNdecisionData)
ANNanalyseOneHistogram('Amoxicillin','Salmonella enterica',newHistosANNdecisionData)
ANNanalyseOneHistogram('Amoxicillin','Streptococcus pneumoniae',newHistosANNdecisionData)
ANNanalyseOneHistogram('Amikacin','Haemophilus influenzae',newHistosANNdecisionData)
ANNanalyseOneHistogram('Amikacin','Mycobacterium tuberculosis',newHistosANNdecisionData)
ANNanalyseOneHistogram('Amikacin','Staphylococcus aureus',newHistosANNdecisionData)

%fungal with no CLSI bp:
ANNanalyseOneHistogram('AmphotecerinB','Candida albicans EUCAST',newHistosANNdecisionData)
ANNanalyseOneHistogram('AmphotecerinB','Candida glabrata EUCAST',newHistosANNdecisionData)
ANNanalyseOneHistogram('Fluconazole','Candida glabrata',newHistosANNdecisionData)


%%

% Put ANN decisions on the spectrally ordered histogram plots
illustrateNewHistoDecisions(newHistosANNdecisionData);

%%

net = fitnet([15 15]);
%net = feedforwardnet([30 30 30 30]);
traintype = 0;
Nbumps = 3;
NNdata = makeSIRTrainingTablesFromNbumps(consolidatedSIRData,Nbumps,traintype);
disp(NNdata);
net = train(net,NNdata.inputTable',NNdata.targetTable');
%net = newrb(NNdata.inputTable',NNdata.targetTable');
