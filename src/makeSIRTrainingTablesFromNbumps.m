function trainingOutput = makeSIRTrainingTablesFromNbumps(SIRData,Nbumps,traintype,ignoreList)

    if nargin < 3
        traintype = 1;
    end
    if nargin < 4
        ignoreList = [];
    end

    close all

    data = extractSIRdataColumns(SIRData);
    drugs = data.drugs;
    bugs = data.bugs;
    Nclusters = data.Nclusters;
    GMMconverged = data.GMMconverged;
    GMMRsquared = data.GMMRsquared;
    GMMSmean = data.GMMSmean;
    GMMSsigma = data.GMMSsigma;
    GMMImean = data.GMMImean;
    GMMIsigma = data.GMMIsigma;
    GMMRmean = data.GMMRmean;
    GMMRsigma = data.GMMRsigma;
    SRboundary = data.SRboundary;
    SIboundary = data.SIboundary;
    IRboundary = data.IRboundary;
    Tecoff = data.Tecoff;
    CLSIsticS = data.CLSIsticS;
    CLSIsticR = data.CLSIsticR;
    Xecoff = data.Xecoff;
    XecoffCorrelation = data.XecoffCorrelation;
    GMMweights = data.GMMweights;

    clear data

    targetTable = [CLSIsticS , CLSIsticR];
    targetVariables = {'CLSISTICS','CLSISTICR'};

    %The N-bump stuff:
    F = (Nclusters == Nbumps);
    F(ignoreList) = 0;
    F = F & GMMconverged;
    F = F & ~isnan(targetTable(:,1));
    F = F & ~isnan(targetTable(:,2));

    switch traintype
        case 0
            trainingOutput.name = 'GMM data with as much info as possible';
            switch Nbumps
                case 3
                    inputTable = ...
                    [GMMSmean(F),sqrt(GMMSsigma(F)),GMMImean(F),sqrt(GMMIsigma(F)),...
                    GMMRmean(F),sqrt(GMMRsigma(F)),...
                    GMMweights(F,1),GMMweights(F,2),GMMweights(F,3),...
                    SIboundary(F),IRboundary(F),SRboundary(F)];
        
                    inputVariables = ...
                        {'Smean','Sstd','Imean','Istd','Rmean','Rstd',...
                        'Sweight','Iweight','Rweight',...
                        'SIboundary','IRboundary','SRboundary'};
                case 2
                    inputTable = ...
                    [GMMSmean(F),sqrt(GMMSsigma(F)),...
                    GMMRmean(F),sqrt(GMMRsigma(F)),...
                    GMMweights(F,1),GMMweights(F,3),SRboundary(F)];
        
                    inputVariables = ...
                        {'Smean','Sstd','Rmean','Rstd',...
                        'Sweight','Rweight','SRboundary'};
                case 1
                    inputTable = [GMMSmean(F),sqrt(GMMSsigma(F))];        
                    inputVariables = {'Smean','Sstd'};
            end
        case 1
            trainingOutput.name = 'all GMM data but without bump weights';
            switch Nbumps
                case 3
                    inputTable = ...
                    [GMMSmean(F),sqrt(GMMSsigma(F)),GMMImean(F),sqrt(GMMIsigma(F)),...
                    GMMRmean(F),sqrt(GMMRsigma(F)),...
                    SIboundary(F),IRboundary(F),SRboundary(F)];
        
                    inputVariables = ...
                        {'Smean','Sstd','Imean','Istd','Rmean','Rstd',...
                        'SIboundary','IRboundary','SRboundary'};
                case 2
                    inputTable = ...
                    [GMMSmean(F),sqrt(GMMSsigma(F)),...
                    GMMRmean(F),sqrt(GMMRsigma(F)),SRboundary(F)];
        
                    inputVariables = ...
                        {'Smean','Sstd','Rmean','Rstd','SRboundary'};
                case 1
                    inputTable = [GMMSmean(F),sqrt(GMMSsigma(F))];        
                    inputVariables = {'Smean','Sstd'};
            end
        case 2
            trainingOutput.name = 'GMM data with only bump locations and variances';

            switch Nbumps
                case 3
                    inputTable = ...
                    [GMMSmean(F),sqrt(GMMSsigma(F)),...
                     GMMImean(F),sqrt(GMMIsigma(F)),...
                     GMMRmean(F),sqrt(GMMRsigma(F))];
                    inputVariables = {'Smean','Sstd','Imean','Istd','Rmean','Rstd'};
                case 2
                    inputTable = ...
                    [GMMSmean(F),sqrt(GMMSsigma(F)),...
                     GMMRmean(F),sqrt(GMMRsigma(F))];
                    inputVariables = {'Smean','Sstd','Rmean','Rstd'};
                case 1
                    inputTable = [GMMSmean(F),sqrt(GMMSsigma(F))];
                    inputVariables = {'Smean','Sstd'};
            end
        case 3
            trainingOutput.name = 'GMM data with only bump boundaries';

            switch Nbumps
                case 3
                    inputTable = [SIboundary(F),IRboundary(F),SRboundary(F)];
                    inputVariables = {'SIboundary','IRboundary','SRboundary'};
                case 2
                    inputTable = [SRboundary(F)];
                    inputVariables = {'SRboundary'};
                case 1
                    inputTable = [];
                    inputVariables = {};
            end
    end

    targetTable = targetTable(F,:);

    if any(isnan(inputTable(:)))
        error('There is NaN in the ANN input data')
    end
    if any(isnan(targetTable(:)))
        error('There is NaN in the ANN target data')
    end
    
    trainingOutput.inputTable = inputTable;
    trainingOutput.inputVariables = inputVariables;
    trainingOutput.targetTable = targetTable;
    trainingOutput.targetVariables = targetVariables;
    
    disp('----------------------------------------------------')
    disp(trainingOutput.name);
    disp('----------------------------------------------------')
    

end