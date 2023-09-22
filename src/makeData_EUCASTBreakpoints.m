%This file is no longer used: REB 19 Feb 2023
 
clear all
close all
clc

%%

EUCASTBreakpoints = importEUCASTBreakpoints('./EUCAST-CLSI-Matching/old-matching-Docs/EUCAST Breakpoints Tidy.csv');
save('./mat/EUCASTbreakpoints.mat','EUCASTBreakpoints')

%%

load('./mat/extractedEUCAST-SIR-Data.mat')

%%

consolidatedCLSIEUCASTMatches = findAllCLSIEUCASTMatches(consolidatedSIRData,EUCASTBreakpoints);
save('./mat/EUCASTCLSImatches.mat','consolidatedCLSIEUCASTMatches')

