%% makePaperDataFiles

clc
clear all
close all

%%

disp('Working on GM modelling....')
makeData_GMMCSVToMAT;

% uses files   : CLSI/breakpoints.csv
%              : EUCASThistograms/*.csv
% makes files  : csv/',thisDrugSp,'.csv' on line 344
%              : csv/_entireDataset.csv
%              : csv/_consolidatedDataset.csv
%              : mat/extractedEUCAST-SIR-Data.mat

%% makeData_EUCASTToMat is run by makeData_GMMCSVToMAT (see last line)

% makeData_EUCASTToMat;

% uses files   : EUCASThistograms/*.csv
% makes files  : csv/_all-EUCAST-Histogtrams.csv
%              : mat/all-EUCAST-Histograms.mat

%%

%makeData_EUCASTBreakpoints;

% uses files   : EUCAST-CLSI-Matching/old-matching-Docs/EUCAST Breakpoints Tidy.csv
%              : mat/extractedEUCAST-SIR-Data.mat
% makes files  : mat/EUCASTbreakpoints.mat
%              : mat/EUCASTCLSImatches.mat

% REB 19 Dec 2022:
% these mat files have been superceded by the next section,
% do not use them anymore. For instance
% consolidatedCLSIEUCASTMatches in EUCASTCLSIMatches.mat
% is out of date.
% It can be used via plotCLSIEUCASTcorrelation(consolidatedCLSIEUCASTMatches);
% but this is out of date usage.

%%

disp('Working on breakpoints....')
makeData_EUCASTBreakpoints2;

% uses files   : EUCAST-CLSI-Matching/df_matched_bp.xlsx
% makes files  : mat/bothBreakpointSets.mat

%%

disp('Incorporating breakpoints....')
incorporateEUCASTbreakpointsIntoSIR;

% uses ./mat/extractedEUCAST-SIR-Data.mat
% modifies variables consolidatedSIRData and entireSIRData
% to incorporate EUCAST breakpoints (but only 232 of them)

%%

clc
clear all
close all

load('./mat/all-EUCAST-Histograms.mat')
load('./mat/extractedEUCAST-SIR-Data.mat')
load('./neuralNets/ebestNets.mat')
load('./neuralNets/cbestNets.mat')

% this saves ./mat/ANNdecisions.mat with the variable ANNdecisionData
% ----------------------------------------------------------------------------------------
% It uses the ANNs in bestNets to replicate CLSI decisions using EUCAST
% histograms, but only for those histograms where CLSI breakpoints exist.
% Those decisions are then held in ANNdecisionData.table:

makeANNcomparisonData(allEUCASTHistograms,consolidatedSIRData,ebestNets,0);%eucast
makeANNcomparisonData(allEUCASTHistograms,consolidatedSIRData,cbestNets,1);%clsi

%%

load('./mat/cANNdecisions.mat')
load('./mat/eANNdecisions.mat')

% this saves ./mat/allnewHistosANNdecisions.mat with the variable newHistosANNdecisionData
% ----------------------------------------------------------------------------------------
% ANN breakpoint decisions are now made for histograms with more than P.numberOfPatientsLimit in
% globalParameterValues.m and those decisions are found in the newHistosANNdecisionData.table.
% Post-decisions linear correction to those using prior validated ANN
% decisions is currently not implemented in the following code
% because it seems that high ANN decisions can be pushed far too high:

determineAllHistogramBreakpoints(allEUCASTHistograms,cANNdecisionData,cbestNets,1)%clsi
determineAllHistogramBreakpoints(allEUCASTHistograms,eANNdecisionData,ebestNets,0)%eucast

load('./mat/c_allnewHistosANNdecisions.mat')
load('./mat/e_allnewHistosANNdecisions.mat')

%%

% this merges breakpoints into the table of new ANN decisions and it
% re-saves ./mat/allnewHistosANNdecisions.mat after updating the variable newHistosANNdecisionData

close all
mergeValidatePriorNewANNdecisions(cANNdecisionData,CnewHistosANNdecisionData,1);%clsi

close all
mergeValidatePriorNewANNdecisions(eANNdecisionData,EnewHistosANNdecisionData,0)%;eucast

%Do not run any other codes at this point .... variables have not been updated
%better to clear and re-load all variables...

clear all

disp('--------')
disp('Finished')
disp('--------')

