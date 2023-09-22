function output = trainNNonAllHistogramData(type,times,CEflag)

    if nargin < 3
        CEflag = 1;
    end
    if nargin < 1
        type = 1;
    end
    if nargin < 2
        times = 1;
    end
    
    % collate data for training a 2-output (S+R breakpoints) neural network
    % where inputs are all the EUCAST histograms. Breakpoints are currently
    % CLSI's.
    
    close all
    
    %%
    
    load('./mat/all-EUCAST-Histograms.mat')
    load('./mat/extractedEUCAST-SIR-Data.mat')
    
    entirePDFs = {entirePDFs{3:end}};
    
    %%
    
    [N,M] = size(allEUCASTHistograms);
    
    inputTable = zeros(N,length(lDosages));
    targetTable = zeros(N,2);
    PApairs = cell(N,2);
    
    output.inputVariableNames = lDosagesString(3:end);
    
    if CEflag
        J = 19;
        K = 20;
        output.targetVariableNames = {'CLSI STIC S','CLSI STIC R'};
    else
        J = 23;
        K = 24;
        output.targetVariableNames = {'EUCAST S','EUCAST R'};
    end
    
    for n = 1:N
        iT = cell2mat(allEUCASTHistograms(n,3:M));
        %iT = iT / sum(iT);
        inputTable(n,:) = iT;
        targetTable(n,:) = cell2mat(entireSIRData(n+2,J:K));
        PApairs(n,1) = entireSIRData(n+2,1);
        PApairs(n,2) = entireSIRData(n+2,2);
    
        checkPApair = strcmpi(allEUCASTHistograms{n,1},entireSIRData{n+2,1}) & ...
            strcmpi(allEUCASTHistograms{n,2},entireSIRData{n+2,2});
        if ~checkPApair
            disp('---------')
            disp(num2str(n));
            disp('---------')
            disp(['|',allEUCASTHistograms{n,1},'|=?|',entireSIRData{n+2,1},'|'])
            disp(['|',allEUCASTHistograms{n,2},'|=?|',entireSIRData{n+2,2},'|'])
            error('A problem needs resolving!')
        end
    end
    
    %%
    
    %At this point lots of target variables are NaN, so remove them:
    INT = isnan(targetTable);
    INT = INT(:,1) | INT(:,2);
    INT = ~INT;

    if ~CEflag
        noIssue = (targetTable(:,1) > -9);
        INT = INT & noIssue;
    end
    
    cleanInputTable = inputTable(INT,:);
    cleanTargetTable = targetTable(INT,:);
    cleanPApairs = PApairs(INT,:);
    cleanPDFs = {entirePDFs{INT}};
    
    [N,M] = size(cleanInputTable);
    smoothcleanInputTable = cleanInputTable;

    disp(['ANN training based on ',num2str(N),' MIC histograms'])

    for n = 1:N
        cIT = cleanInputTable(n,:);
        cIT = cIT / sum(cIT);
        cleanInputTable(n,:) = cIT;
    
        smoothcIT = cleanPDFs{n}(lDosages');
        smoothcleanInputTable(n,:) = smoothcIT;
    
        %plot(smoothcIT)
        %hold on
        %plot(cIT,'ok')
        %pause
        %close(1)
    end
    
    %%
    
    %nftool;
    
    %Build a NN based on inputs and outputs:
    
    switch type
        case 1
            trainVector = [30 30 30 30];
            net = feedforwardnet(trainVector);
            netfun = @feedforwardnet;
            NNname = 'feedforward net';
        case 2
            %Radial basis function ANN: is this just interpolation?
            disp('Radial basis net specified, so only training once ...')
            net = newrb(smoothcleanInputTable',cleanTargetTable');
            netfun = @newrb;
            NNname = 'radial basis net';
            netpredictTarget = net(smoothcleanInputTable');
            corrPerformances = regression(cleanTargetTable',netpredictTarget,'one');
            rmsPerformances = perform(net,netpredictTarget,cleanTargetTable');
            trainVector = NaN;
            trainVectorNbrd = NaN;
            times = 1;
        case 3
            trainVector = [50 50];
            net = fitnet(trainVector);
            netfun = @fitnet;
            NNname = 'matlab fitnet';
        case 4
            trainVector = [20 20 20];
            net = feedforwardnet(trainVector);
            netfun = @feedforwardnet;
            NNname = 'feedforward net';
        case 5
            trainVector = [30 30];
            net = feedforwardnet(trainVector);
            netfun = @feedforwardnet;
            NNname = 'feedforward net';
        case 6
            trainVector = [40 40];
            net = cascadeforwardnet(trainVector);
            netfun = @cascadeforwardnet;
            NNname = 'cascadeforward net';
    end
    
    if ~(type == 2)
        trainVectorNbrd = vectorNeighbour(trainVector,times);
        [ntv,~] = size(trainVectorNbrd);
        %reduce times if need be:
        times = ntv;
        %for j = 1:length(net.layers)
            %net.layers{j}.transferFcn = 'tansig';
            %net.layers{j}.transferFcn = 'logsig';
            %net.layers{j}.transferFcn = 'elliotsig';
        %end
        
        %initialise variables:
        nets = cell(times,1);
        corrPerformances = zeros(1,times);
        rmsPerformances = zeros(1,times);
        netpredictTarget = cell(times,1);
    
        for n = 1:times
            nets{n} = net;
        end
        disp('--------------------')
        disp(['Training a ',NNname])
        disp('--------------------')
        for n = 1:times
            disp(['NN training ',num2str(100*n/times,3),'% complete ...']);
            disp(['Neuron structure parameters : ',num2str(trainVectorNbrd(n,:))]);
            net = netfun(trainVectorNbrd(n,:));
            if times > 1
                net.trainParam.showWindow = 0;
            end
            net = train(net,smoothcleanInputTable',cleanTargetTable','UseParallel','yes');
            nets{n} = net;
            netpredictTarget{n} = net(smoothcleanInputTable');
            corrPerformances(n) = regression(cleanTargetTable',netpredictTarget{n},'one');
            rmsPerformances(n) = perform(net,netpredictTarget{n},cleanTargetTable');
            disp(['Performance where (1,0) is best : ',num2str([corrPerformances(n),rmsPerformances(n)])]);
        end
    
        output.N = times;
        output.nets = nets;
    
    else
        
        close all
        output.N = 1;
        output.nets = net;
    
    end
    
    output.NNstructures = trainVectorNbrd;
    output.netfun = netfun;
    output.corrPerformances = corrPerformances;
    output.rmsPerformances = rmsPerformances;
    output.name = NNname;
    output.netpredictTarget = netpredictTarget;
    
    [~,jc] = max(corrPerformances);
    
    if times > 1
        %re-train and visualise the best NN:
        net = nets{jc};
        net.trainParam.showWindow = 1;
        net = train(net,smoothcleanInputTable',cleanTargetTable','useParallel','yes');
    
        netpredictTarget{jc} = net(smoothcleanInputTable');
        output.netpredictTarget{jc} = netpredictTarget{jc};
    
        output.corrPerformances(jc) = regression(cleanTargetTable',netpredictTarget{jc},'one');
        output.rmsPerformances(jc) = perform(net,netpredictTarget{jc},cleanTargetTable');
        output.nets{jc} = net;
    end
    
    output.targetData = cleanTargetTable';
    output.trainingData = smoothcleanInputTable';
    
    %disp('Plotting correlations for best performing network')
    %if times > 1
    %    plotregression(cleanTargetTable',netpredictTarget{jc})
    %end
    %view(net)

    disp(['training data ratio ',num2str(net.divideParam.trainRatio)]);
    disp(['test data ratio ',num2str(net.divideParam.testRatio)]);
    disp(['validation data ratio ',num2str(net.divideParam.valRatio)]);
    
end
