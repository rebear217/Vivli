function ANNdecisionData = makeANNcomparisonData(allEUCASTHistograms,consolidatedSIRData,bestNetsEnsemble,CEflag)

    if nargin < 4
        CEflag = 1;
    end

    close all

    data = extractSIRdataColumns(consolidatedSIRData);
    lDosages = data.lDosages;
    %lDosagesString = data.lDosagesString;

    if CEflag
        breakpointsS = data.CLSIsticS;
        breakpointsR = data.CLSIsticR;
        bpSname = 'CLSIbpS';
        bpRname = 'CLSIbpR';
    else
        breakpointsS = data.eucastBPS;
        breakpointsR = data.eucastBPR;
        bpSname = 'EUCASTbpS';
        bpRname = 'EUCASTbpR';
    end

    drugs = data.drugs;
    bugs = data.bugs;

    clear data

    [n,~] = size(consolidatedSIRData);
    histoPos = 1;
    histoPositions = NaN(1,n);
    histoPositions(1) = -1;
    histoPositions(2) = -1;

    for j = 3:n
        consString = [consolidatedSIRData{j,1},consolidatedSIRData{j,2}];
        histoString = [allEUCASTHistograms{histoPos,1},allEUCASTHistograms{histoPos,2}];
        while ~strcmpi(consString,histoString)
            histoPos = histoPos + 1;
            histoString = [allEUCASTHistograms{histoPos,1},allEUCASTHistograms{histoPos,2}];
        end
        histoPositions(j) = histoPos;
    end
    if any(isnan(histoPositions))
        error('Cannot match all consolidated SIR data drug-bug names into the all EUCAST histograms variable')
    end

    %histoPositions is for referencing into allEUCASTHistograms
    
    %usePositions is for referencing into consolidatedSIRData
    %for where breakpoints are not NaN.

    %So histoPositions(usePositions) says where the histograms are
    %that have defined breakpoints

    histoPositions = histoPositions(3:end);
    allHistograms = cell2mat(allEUCASTHistograms(:,3:end));
    numberOfPatients = sum(allHistograms,2);
    HR = repmat(numberOfPatients,1,19);

    allHistograms = allHistograms ./ HR;

    usePositions = find(~isnan(breakpointsS));
    checkusePositions = find(~isnan(breakpointsR));
    if any(usePositions & ~checkusePositions)
        warning('breakpoint S and R have different NaNs')
    else
        disp('check: breakpoint S and R have the same NaNs (which is good)')
    end
    
    N = length(usePositions);
    disp(['Estimating breakpoints for ',num2str(N),' strains using the ANN ensemble ...'])

    Drugs = cell(N,1);
    Bugs = cell(N,1);

    bpS = breakpointsS(usePositions);
    bpR = breakpointsR(usePositions);

    ANNbpS = zeros(size(bpS));
    ANNbpR = zeros(size(bpR));
    ANNbpSIQR = zeros(length(bpS),2);
    ANNbpRIQR = zeros(length(bpR),2);    
    ANNbpSDecisions = zeros(N,bestNetsEnsemble.N);
    ANNbpRDecisions = zeros(N,bestNetsEnsemble.N);

    for j = 1:N

        Drugs{j} = drugs{usePositions(j)};
        Bugs{j} = bugs{usePositions(j)};

        if mod(j,10) == 0
            disp([num2str(100*j/N,3),'% complete'])
        end
        
        J = histoPositions(usePositions(j));
        HJ = allHistograms(J,:);
        ANNbreakpointDecisions = ANNdecision(HJ',bestNetsEnsemble,[],CEflag);

        ANNbpS(j) = ANNbreakpointDecisions.bpLeftFinal;
        ANNbpR(j) = ANNbreakpointDecisions.bpRightFinal;
        ANNbpRIQR(j,:) = ANNbreakpointDecisions.bfRightIQR;
        ANNbpSIQR(j,:) = ANNbreakpointDecisions.bfLeftIQR;
        ANNbpSDecisions(j,:) = ANNbreakpointDecisions.bpLefts;
        ANNbpRDecisions(j,:) = ANNbreakpointDecisions.bpRights;    

    end

    mdlS = fitlm(bpS,ANNbpS);
    aS = mdlS.Coefficients.Estimate(1);
    bS = mdlS.Coefficients.Estimate(2);
    correctS = @(x)((x-aS)/bS);

    mdlR = fitlm(bpR,ANNbpR);
    aR = mdlR.Coefficients.Estimate(1);
    bR = mdlR.Coefficients.Estimate(2);
    correctR = @(x)((x-aR)/bR);

    ANNdecisionData.table = table(Drugs,Bugs,bpS,bpR,...
        ANNbpS,ANNbpR,...
        ANNbpSIQR(:,1),ANNbpSIQR(:,2),...
        ANNbpRIQR(:,1),ANNbpRIQR(:,2),...
        'VariableNames',{'drugs','bugs',bpSname,bpRname,'ANNbpS','ANNbpR','ANNbpSIQRlow','ANNbpSIQRhigh','ANNbpRIQRlow','ANNbpRIQRhigh'});
    
    ANNdecisionData.correctedTable = table(Drugs,Bugs,bpS,bpR,...
        correctS(ANNbpS),correctR(ANNbpR),...
        correctS(ANNbpSIQR(:,1)),correctS(ANNbpSIQR(:,2)),...
        correctR(ANNbpRIQR(:,1)),correctR(ANNbpRIQR(:,2)),...
        'VariableNames',{'drugs','bugs',bpSname,bpRname,'ANNbpS','ANNbpR','ANNbpSIQRlow','ANNbpSIQRhigh','ANNbpRIQRlow','ANNbpRIQRhigh'});

    ANNdecisionData.ensembleSDecisions = ANNbpSDecisions;
    ANNdecisionData.ensembleRDecisions = ANNbpRDecisions;    
    ANNdecisionData.note = 'In ensemble{S/R}Decisions, columns are the ANN decisions and rows are drug-bug pairs commensurate with the table ordering';
    ANNdecisionData.histograms = allHistograms(histoPositions(usePositions),:);
    ANNdecisionData.numberOfPatients = numberOfPatients(histoPositions(usePositions));

    if CEflag
        cANNdecisionData = ANNdecisionData;
        save('./mat/cANNdecisions.mat','cANNdecisionData')
    else
        eANNdecisionData = ANNdecisionData;
        save('./mat/eANNdecisions.mat','eANNdecisionData')
    end

end