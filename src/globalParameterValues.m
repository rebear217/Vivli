function P = globalParameterValues()
    
    %Use when making synthetic histograms:
    P.noiseLevel = 0.0;

    %this is used in linear modelling aspects of analyseSIR2bump
    P.ignoreDiscrepancy = 2;
    P.Replicates = 20;

    %max number of GMM models to use:
    P.Ngmmodels = 3;

    %normal tail parameter of where to place the ECOFF:
    P.ECOFFpValue = 0.99;

    P.lDosages = [-8.9658 -7.9658 -6.9658 -5.9658 -5.0589 -4.0589 -3 -2 -1 0 1 2 3 4 5 6 7 8 9];
    P.lDosagesString = {'-8.9658','-7.9658','-6.9658','-5.9658','-5.0589','-4.0589','-3','-2','-1','0','1','2','3','4','5','6','7','8','9'};
    
    % used in ANN ensemble final decision making
    P.minQRtile = 25;
    P.maxQRtile = 75;
    P.outlierValue = 25;
    
    %minimal limit of EUCAST histo size to determine new cutoffs with ANNs:
    P.numberOfPatientsLimit = 100;
    %This value leads to 482 valid drug-bug pairs with CLSI breakpoints in newHistosANNdecisionData.table
    %whereas there are 590 such entries in ANNdecisionData.table
    
    %Linear correction post-processing made things worse, turn it off:
    P.performLinearCorrectionPostProcessing = 0;
    %This was implemented to cope with the apparent systematic error that
    %low MICs were over-estimated and high MICs were under-estimated by the
    %ANNs ensemle. But correcting this makes other metrics worse.
    %Obviously!

end