function ANNbreakpointDecisions = ANNdecision(inputDataInstance,bestNetsEnsemble,priorANNdecisionData,CEflag)

    if nargin < 3
        priorANNdecisionData = [];
    end

    if nargin < 4
        CEflag = 1;
    end
    
    P = globalParameterValues();
    breakpoint = zeros(2,bestNetsEnsemble.N);

    G = gcp;
    if ~G.Connected
        disp('Starting parallel toolbox for ANN ensemble evaluation')
    end

    parfor k = 1:bestNetsEnsemble.N
        breakpoint(:,k) = bestNetsEnsemble.nets{k}(inputDataInstance);
    end

    if ~isempty(priorANNdecisionData)

        %once we have made and tested prior ANN decisions, we can compare the ensemble
        %for correctness and shift bad decisions in better directions with a correction:

        if CEflag
            mdlS = fitlm(priorANNdecisionData.table.CLSIbpS,priorANNdecisionData.table.ANNbpS);
        else
            mdlS = fitlm(priorANNdecisionData.table.EUCASTbpS,priorANNdecisionData.table.ANNbpS);
        end
        aS = mdlS.Coefficients.Estimate(1);
        bS = mdlS.Coefficients.Estimate(2);
        correctS = @(x)((x-aS)/bS);
    
        if CEflag
            mdlR = fitlm(priorANNdecisionData.table.CLSIbpR,priorANNdecisionData.table.ANNbpR);
        else
            mdlR = fitlm(priorANNdecisionData.table.EUCASTbpR,priorANNdecisionData.table.ANNbpR);
        end
        aR = mdlR.Coefficients.Estimate(1);
        bR = mdlR.Coefficients.Estimate(2);
        correctR = @(x)((x-aR)/bR);

    else
        %if we have no prior decisions to test, no corrections can be made, so "correct" using
        %identity functions:
        correctR = @(x)x;
        correctS = @(x)x;
    end

    ANNbreakpointDecisions.bpLefts = correctS(breakpoint(1,:));
    ANNbreakpointDecisions.bpRights = correctR(breakpoint(2,:));

    decisionIQR = @(x)prctile(x,[P.minQRtile P.maxQRtile]);
    decision = @(x)trimmean(x,P.outlierValue);
    %decision = @(x)median(x);

    ANNbreakpointDecisions.bpLeftFinal = decision(ANNbreakpointDecisions.bpLefts);
    ANNbreakpointDecisions.bpRightFinal = decision(ANNbreakpointDecisions.bpRights);

    ANNbreakpointDecisions.bfLeftIQR = decisionIQR(ANNbreakpointDecisions.bpLefts);
    ANNbreakpointDecisions.bfRightIQR = decisionIQR(ANNbreakpointDecisions.bpRights);

end
