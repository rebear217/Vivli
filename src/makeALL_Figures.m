%% makePaperfigures

cd('/Users/RobertBeardmore/Library/CloudStorage/Dropbox/EUCAST-ECOFF')

clc
clear all
close all

%%

load('./mat/all-EUCAST-Histograms.mat')
load('./mat/extractedEUCAST-SIR-Data.mat')
load('./mat/EUCASTCLSImatches.mat');
load('./neuralNets/ebestNets.mat')
load('./neuralNets/cbestNets.mat')
load('./mat/bothBreakpointSets.mat');
load('./mat/e_allnewHistosANNdecisions.mat')
load('./mat/c_allnewHistosANNdecisions.mat')
load('./mat/eANNdecisions.mat')
load('./mat/cANNdecisions.mat')
load('./mat/weightVector.mat')
load('./mat/orderedUSEUdifferencesTable.mat')
rawCLSIbreakpoints = importCLSIbreakpoint('./CLSI/breakpoints.csv');
showbreakpointBreakdown;
load('./mat/AMR-R.mat')

%%

% this checks for errors between SIR variables and the bothBreakpointSets variables, updating in case of NaNs:
% It has not reported any errors, but it has reported NaNs that can be udpated.

[bothBreakpointSets,consolidatedSIRData] = validateSIRdata(bothBreakpointSets,consolidatedSIRData);
[~,entireSIRData] = validateSIRdata(bothBreakpointSets,entireSIRData);

%% Figure 1

close all
plotCLSIEUCASTcorrelation2(bothBreakpointSets);
figure(1)
export_fig('./figures/1.pdf')
figure(3)
export_fig('./figures/3.pdf')
figure(4)
export_fig('./figures/4.pdf')

%%

close all

%CLSI breakpoints:
spectralEmbedding(allEUCASTHistograms,consolidatedSIRData,1,weight);

figure(2)
export_fig('./figures/5.pdf')
figure(4)
export_fig('./figures/6.pdf')
print('-r600','-dpng','./figures/6.png')
figure(5)
export_fig('./figures/7.pdf')
figure(3)
export_fig('./figures/8.pdf')
figure(19)
export_fig('./figures/9.pdf')
figure(75)
export_fig('./figures/10.pdf')
figure(79)
export_fig('./figures/10_.pdf')

%%

%EUCAST breakpoints:
spectralEmbedding(allEUCASTHistograms,consolidatedSIRData,0,weight);

figure(2)
export_fig('./figures/5_e.pdf')
figure(4)
export_fig('./figures/6_e.pdf')
print('-r600','-dpng','./figures/6_e.png')
figure(5)
export_fig('./figures/7_e.pdf')
figure(3)
export_fig('./figures/8_e.pdf')
figure(42)
export_fig('./figures/9_e.pdf')
figure(62)
export_fig('./figures/10_e.pdf')


%% Figure 2

close all
ECOFFexamples();
figure(1)
export_fig('./figures/11.pdf')
figure(2)
export_fig('./figures/12.pdf')

close all
analyseOneHistogram(allEUCASTHistograms,'Cefotaxime','Klebsiella aerogenes')
figure(1)
export_fig('./figures/13.pdf')
figure(2)
export_fig('./figures/14.pdf')

%%

close all
summariseSIRdata(consolidatedSIRData,1);

figure(1)
set(1,'pos',[55 1 2506 696]);
xlim([1 439])
export_fig('./figures/16.pdf')

figure(2)
xlim([149        200])
export_fig('./figures/17.pdf')

figure(5)
export_fig('./figures/18.pdf')

figure(6)
export_fig('./figures/18_.pdf')

%%

close all
CEflag = 0;
summariseSIRdata(consolidatedSIRData,CEflag);

figure(1)
set(1,'pos',[55 1 2506 696]);
xlim([1 348])
export_fig('./figures/16e.pdf')

figure(2)
xlim([290       348])
export_fig('./figures/17e.pdf')

figure(6)
export_fig('./figures/18_e.pdf')

%% Figure 4

close all
CFflag = 1;
illustrateNewHistoDecisions(CnewHistosANNdecisionData,CFflag);%clsi
figure(1)
export_fig('./figures/19.pdf')
figure(2)
export_fig('./figures/20.pdf')
figure(3)
export_fig('./figures/21.pdf')

close all
CFflag = 0;
illustrateNewHistoDecisions(EnewHistosANNdecisionData,CFflag);%eucast
figure(1)
export_fig('./figures/19_e.pdf')
figure(2)
export_fig('./figures/20_e.pdf')
figure(3)
export_fig('./figures/21_e.pdf')


%% Supp Figures...

%% Supp Figure 1

close all

CFflag = 1;
ANNanalyseOneHistogram('Vancomycin','Staphylococcus aureus',cANNdecisionData,CFflag,1)
export_fig('./figures/22.pdf')
ANNanalyseOneHistogram('Ertapenem','Escherichia coli',cANNdecisionData,CFflag,3)
export_fig('./figures/24.pdf')
ANNanalyseOneHistogram('Fluconazole','Candida glabrata',CnewHistosANNdecisionData,CFflag,5)
export_fig('./figures/25.pdf')

close all

CFflag = 0;
ANNanalyseOneHistogram('Vancomycin','Staphylococcus aureus',eANNdecisionData,CFflag,1)
export_fig('./figures/22_e.pdf')
ANNanalyseOneHistogram('Ertapenem','Escherichia coli',eANNdecisionData,CFflag,3)
export_fig('./figures/24_e.pdf')
ANNanalyseOneHistogram('Fluconazole','Candida glabrata',EnewHistosANNdecisionData,CFflag,5)
export_fig('./figures/25_e.pdf')

%% Supp Figure 2

close all

analyseSIR1bump(consolidatedSIRData);
figure(4)
export_fig('./figures/26.pdf')
close all

analyseSIR2bump(consolidatedSIRData,1);
figure(7)
export_fig('./figures/27.pdf')
close all

analyseSIR3bump(consolidatedSIRData,1);
figure(4)
export_fig('./figures/28.pdf')
close all

%%

close all
CEflag = 1;
testBestNNs(cbestNets,CEflag,lDosages,2);
figure(1)
export_fig('./figures/33.pdf')
figure(2)
export_fig('./figures/34.pdf')

close all
CEflag = 0;
testBestNNs(ebestNets,CEflag,lDosages,2);
figure(1)
export_fig('./figures/33_e.pdf')
figure(2)
export_fig('./figures/34_e.pdf')

%%

close all
gmmFitsHistogramPlot(entireSIRData)
export_fig('./figures/35.pdf')

%%

close all
CLSIenterobacteralesCount(bothBreakpointSets,rawCLSIbreakpoints)
figure(1)
export_fig('./figures/36.pdf')

%%

close all
compareANNEs(CnewHistosANNdecisionData,EnewHistosANNdecisionData);
figure(1)
export_fig('./figures/37.pdf')
figure(3)
export_fig('./figures/38.pdf')

%%

close all
orderedUSEUdifferencesTable = plotBreakpointDifferences(bothBreakpointSets);
save('./mat/orderedUSEUdifferencesTable.mat','orderedUSEUdifferencesTable')
writetable(orderedUSEUdifferencesTable,'./csv/orderedUSEUdifferencesTable.csv')  
figure(2)
export_fig('./figures/39.pdf')

%%

% Do this first to get the weight vector:

close all
weight = weightedRegressionAnalysis(allEUCASTHistograms,consolidatedSIRData);
save('./mat/weightVector.mat','weight')

figure(1)
export_fig('./figures/40.pdf')
figure(2)
export_fig('./figures/41.pdf')
figure(3)
export_fig('./figures/42.pdf')

%%

close all
DeOptimiseQ(weight,cbestNets,ebestNets);

figure(1)
export_fig('./figures/43.pdf')
figure(2)
export_fig('./figures/44.pdf')
figure(3)
export_fig('./figures/45.pdf')
figure(4)
export_fig('./figures/46.pdf')

%%

cd('/Users/RobertBeardmore/Library/CloudStorage/Dropbox/EUCAST-ECOFF/EmilyInnoc')
close all
inocEffect();

%%

close all
CEflag = 0;
ANNanalyseOneHistogram('Levofloxacin','Haemophilus influenzae',eANNdecisionData,CEflag,1)
export_fig('./figures/47.pdf')
CEflag = 1;
ANNanalyseOneHistogram('Levofloxacin','Haemophilus influenzae',cANNdecisionData,CEflag,2)
export_fig('./figures/48.pdf')

%%

close all
processAMRforRlist(AMRforRMatches);
figure(1)
export_fig('./figures/49.pdf')
figure(2)
export_fig('./figures/50.pdf')
figure(3)
export_fig('./figures/51.pdf')

%%

close all
compareANNEecDecisions(eANNdecisionData.table,cANNdecisionData.table)
figure(2)
export_fig('./figures/52_e.pdf')
figure(3)
export_fig('./figures/52_c.pdf')

%%

close all
CEflag = 0;
ANNEsurface(ebestNets,CEflag);
export_fig('./figures/53_e.pdf')

close all
CEflag = 1;
ANNEsurface(cbestNets,CEflag);
export_fig('./figures/53_c.pdf')

%%

cd('/Users/RobertBeardmore/Dropbox/EUCAST-ECOFF')

disp('--------')
disp('Finished')
disp('--------')

close all

