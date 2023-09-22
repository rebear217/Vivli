function determineAllHistogramBreakpoints(allEUCASTHistograms,ANNdecisionData,bestNets,CEflag)

    if nargin < 4
        CEflag = 1;
    end

    if CEflag
        bpSname = 'CLSIbpS';
        bpRname = 'CLSIbpR';
    else
        bpSname = 'EUCASTbpS';
        bpRname = 'EUCASTbpR';
    end

    P = globalParameterValues();

    %[n,~] = size(allEUCASTHistograms);
    %n = n - 2;

    allHistograms = cell2mat(allEUCASTHistograms(:,3:end));
    numberOfPatients = sum(allHistograms,2);
    HR = repmat(numberOfPatients,1,19);
    allHistograms = allHistograms ./ HR;
    
    usePositions = find(numberOfPatients > P.numberOfPatientsLimit);
    n = length(usePositions);
    disp(['Estimating breakpoints for ',num2str(n),' strains using the ANN ensemble ...'])

    ANNbpS = zeros(n,1);
    ANNbpR = zeros(n,1);
    ANNbpSIQR = zeros(n,2);
    ANNbpRIQR = zeros(n,2);    
    ANNbpSDecisions = zeros(n,bestNets.N);
    ANNbpRDecisions = zeros(n,bestNets.N);

    Drugs = cell(n,1);
    Bugs = cell(n,1);

    for J = 1:n

        j = usePositions(J);
        if mod(J,10) == 0
            disp([num2str(100*J/n,3),'% complete'])
        end

        histo = allHistograms(j,:)';

        Drugs{J} = allEUCASTHistograms{j,1};
        Bugs{J} = allEUCASTHistograms{j,2};    
    
        if P.performLinearCorrectionPostProcessing
            %do use prior decision linear corrections:
            ANNbreakpointDecisions = ANNdecision(histo,bestNets,ANNdecisionData,CEflag);
        else
            %do NOT use prior decision linear corrections:
            ANNbreakpointDecisions = ANNdecision(histo,bestNets,[],CEflag);
        end
    
        ANNbpS(J) = ANNbreakpointDecisions.bpLeftFinal;
        ANNbpR(J) = ANNbreakpointDecisions.bpRightFinal;
        ANNbpRIQR(J,:) = ANNbreakpointDecisions.bfRightIQR;
        ANNbpSIQR(J,:) = ANNbreakpointDecisions.bfLeftIQR;
        ANNbpSDecisions(J,:) = ANNbreakpointDecisions.bpLefts;
        ANNbpRDecisions(J,:) = ANNbreakpointDecisions.bpRights;

    end

    bpS = NaN(size(ANNbpS));
    bpR = NaN(size(ANNbpS));

    newHistosANNdecisionData.table = table(Drugs,Bugs,bpS,bpR,ANNbpS,ANNbpR,ANNbpSIQR(:,1),ANNbpSIQR(:,2),ANNbpRIQR(:,1),ANNbpRIQR(:,2),...
        'VariableNames',{'drugs','bugs',bpSname,bpRname,'ANNbpS','ANNbpR','ANNbpSIQRlow','ANNbpSIQRhigh','ANNbpRIQRlow','ANNbpRIQRhigh'});

    newHistosANNdecisionData.ensembleSDecisions = ANNbpSDecisions;
    newHistosANNdecisionData.ensembleRDecisions = ANNbpRDecisions;    
    newHistosANNdecisionData.note = 'In ensemble{S/R}Decisions, columns are the ANN decisions and rows are drug-bug pairs commensurate with the table ordering';
    newHistosANNdecisionData.numberOfPatients = numberOfPatients(usePositions);
    newHistosANNdecisionData.histograms = allHistograms(usePositions,:);

    if CEflag
        CnewHistosANNdecisionData = newHistosANNdecisionData;
        save('./mat/c_allnewHistosANNdecisions.mat','CnewHistosANNdecisionData')
    else
        EnewHistosANNdecisionData = newHistosANNdecisionData;
        save('./mat/e_allnewHistosANNdecisions.mat','EnewHistosANNdecisionData')
    end

end