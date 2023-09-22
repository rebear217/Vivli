function weight = weightedRegressionAnalysis(allEUCASTHistograms,consolidatedSIRData)

    close all

    data = extractSIRdataColumns(consolidatedSIRData);
    lDosages = data.lDosages;
    lDosages = round(lDosages,2);
    bpSclsi = data.CLSIsticS;
    bpRclsi = data.CLSIsticR;
    bpSeucast = data.eucastBPS;
    bpReucast = data.eucastBPR;

    clear data

    weightMatrix = zeros(4,length(lDosages));
    wcount = 1;

    for CEflag = [0 1]

        if CEflag
            breakpointsS = bpSclsi;
            breakpointsR = bpRclsi;
            label = 'CLSI';
        else
            breakpointsS = bpSeucast;
            breakpointsR = bpReucast;
            label = 'EUCAST';
        end
    
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
        HR = repmat(sum(allHistograms,2),1,19);
        allHistograms = allHistograms ./ HR;

        usePositions = find(~isnan(breakpointsS));
        checkusePositions = find(~isnan(breakpointsR));


        if any(usePositions & ~checkusePositions)
            warning('breakpoint S and R have different NaNs')
        else
            disp('check: breakpoint S and R have the same NaNs (which is good)')
        end

        usePositions = intersect(usePositions,find(breakpointsS > -9));
        usePositions = intersect(usePositions,find(breakpointsR > -9));
    
        [paramsS,fitS,linearFitS] = findOptimalHistogramRegression(breakpointsS(usePositions),allHistograms(histoPositions(usePositions),:));
        [paramsR,fitR,linearFitR] = findOptimalHistogramRegression(breakpointsR(usePositions),allHistograms(histoPositions(usePositions),:));

        wcovS = fitS.Coefficients.SE(2:end)./fitS.Coefficients.Estimate(2:end);
        wcovR = fitR.Coefficients.SE(2:end)./fitR.Coefficients.Estimate(2:end);
        fS = (fitS.Coefficients.pValue(2:end) < 0.05);
        errS = wcovS.* paramsS.w;
        fR = (fitR.Coefficients.pValue(2:end) < 0.05);
        errR = wcovR.* paramsR.w;
        
        %now use a further linear fit to map this weighted high-dimensional
        %regression onto better quantitative predictions in a way that
        %won't change the correlation coefficient of the fit:

        XS = breakpointsS(usePositions);
        YS = linearFitS.XweightedFit;
        XR = breakpointsR(usePositions);
        YR = linearFitR.XweightedFit;

        lfitS = fitlm(YS,XS);
        lfitR = fitlm(YR,XR);

        figure(1 + CEflag)
        set(1 + CEflag,'pos',[329         886        1094         800])

        subplot(2,2,1)
        colororder({'k','b'})

        plot(lDosages,paramsS.w,'LineStyle','-','DisplayName','S breakpoint weights');
        hold on
        errorbar(lDosages(fS),paramsS.w(fS),errS(fS),'.','LineStyle','none','DisplayName','S CoV','Color','k','markersize',20);

        xlabel('MIC (log2 \mug/mL)')
        ylabel('explanatory weight (w_j)')
        legend('location','northwest')
        axis tight
        set(gca,'Xtick',lDosages)

        yyaxis right
        plot(lDosages,(fitS.Coefficients.pValue(2:end)),'-b','DisplayName','p-values','LineWidth',0.5);
        ylabel('unitless')

        subplot(2,2,2)
        plot(XS,lfitS.feval(YS),'.k','markersize',24,'DisplayName',['\rho\approx ',num2str(linearFitS.correlation,2)])
        hold on
        xlabel([label,' breakpoints (S, log2 \mug/mL)'])
        ylabel('weighted MIC histogram prediction')
        axis tight
        set(gca,'Xtick',lDosages)
        legend()

        subplot(2,2,3)        
        colororder({'k','b'})
        plot(lDosages,paramsR.w,'LineStyle','-','DisplayName','R breakpoint weights');
        hold on
        errorbar(lDosages(fR),paramsR.w(fR),errR(fR),'.','LineStyle','none','DisplayName','R CoV','Color','k','markersize',20);
        xlabel('MIC (log2 \mug/mL)')
        ylabel('explanatory weight (w_j)')
        legend('location','northwest')
        set(gca,'Xtick',lDosages)
        axis tight

        yyaxis right
        plot(lDosages,(fitR.Coefficients.pValue(2:end)),'-b','DisplayName','p-values','LineWidth',0.5);
        ylabel('unitless')
                
        subplot(2,2,4)
        plot(XR,lfitR.feval(YR),'.k','markersize',24,'DisplayName',['\rho\approx ',num2str(linearFitR.correlation,2)])
        hold on
        xlabel([label,' breakpoints (R, log2 \mug/mL)'])
        ylabel('weighted MIC histogram prediction')
        axis tight
        set(gca,'Xtick',lDosages)
        legend()

        weightMatrix(1 + 2*CEflag,:) = paramsS.w;
        weightMatrix(2 + 2*CEflag,:) = paramsR.w;

    end

    %produce a weight vector by averaging all the weights and then filter
    %noise using moving averaging:

    weight = movmean(mean(weightMatrix,1),4);
    figure(3)
    plot(lDosages,weight,'.-k','markersize',24,'DisplayName','filtered mean (used as w_j)')
    xlabel('MIC (log2 \mug/mL)')
    ylabel('weight (w_j)')
    axis tight
    set(gca,'Xtick',lDosages)    
    hold on
    for j = 1:4
        if j == 4
            plot(lDosages,weightMatrix(j,:),'-','color',[1 1 1]*0.7,'DisplayName','4 weight vectors (CLSI S&R, EUCAST S&R)')
        else
            plot(lDosages,weightMatrix(j,:),'-','color',[1 1 1]*0.7,'HandleVisibility','off')
        end
    end
    legend('location','northwest')
end