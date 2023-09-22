clc
clear all
close all

%%

load('./mat/extractedEUCAST-SIR-Data.mat')
load('./mat/bothBreakpointSets.mat');

%%

consolidatedSIRData = incorporateEUCASTbreakpoints(consolidatedSIRData,bothBreakpointSets);
testCountBreakpointMatchings(consolidatedSIRData,bothBreakpointSets)

%%

entireSIRData = incorporateEUCASTbreakpoints(entireSIRData,bothBreakpointSets);

%%

save('./mat/extractedEUCAST-SIR-Data.mat','consolidatedSIRData','entireSIRData','entirePDFs')
