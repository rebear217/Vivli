function compareANNwithBP(ANNdecisionData,CEflag,figNum,dotColour)

    if nargin < 3
        figNum = 1;
    end
    if nargin < 4
        dotColour = 'r';
    end
    
    P = globalParameterValues();
    if CEflag
        useJ = ~isnan(ANNdecisionData.table.CLSIbpS);
    else
        useJ = ~isnan(ANNdecisionData.table.EUCASTbpS);
        useJ = useJ & (ANNdecisionData.table.EUCASTbpS > -9);
    end

    figure(figNum)
    f = @(p,x)(p(1) + p(2)*x + p(3)*x.^2 + p(4)*x.^3);

    if CEflag
        [rho,pval] = corr(ANNdecisionData.table.CLSIbpS(useJ),ANNdecisionData.table.ANNbpS(useJ),'type','Pearson');
        NLmdlS = fitnlm(ANNdecisionData.table.CLSIbpS(useJ),ANNdecisionData.table.ANNbpS(useJ),f,[0 1 0 2]);
        mdlS = fitlm(ANNdecisionData.table.CLSIbpS(useJ),ANNdecisionData.table.ANNbpS(useJ));
    else
        [rho,pval] = corr(ANNdecisionData.table.EUCASTbpS(useJ),ANNdecisionData.table.ANNbpS(useJ),'type','Pearson');
        NLmdlS = fitnlm(ANNdecisionData.table.EUCASTbpS(useJ),ANNdecisionData.table.ANNbpS(useJ),f,[0 1 0 2]);
        mdlS = fitlm(ANNdecisionData.table.EUCASTbpS(useJ),ANNdecisionData.table.ANNbpS(useJ));
    end

    if CEflag
        plot(ANNdecisionData.table.CLSIbpS(useJ),ANNdecisionData.table.ANNbpS(useJ),'o','color',dotColour,'DisplayName',['raw decisions (r \approx ',num2str(rho,3),')'],'MarkerSize',14)
    else
        plot(ANNdecisionData.table.EUCASTbpS(useJ),ANNdecisionData.table.ANNbpS(useJ),'o','color',dotColour,'DisplayName',['raw decisions (r \approx ',num2str(rho,3),')'],'MarkerSize',14)
    end
    hold on
    try
        if CEflag
            [rhoc,pvalc] = corr(ANNdecisionData.correctedTable.CLSIbpS(useJ),ANNdecisionData.correctedTable.ANNbpS(useJ),'type','Pearson');
            plot(ANNdecisionData.correctedTable.CLSIbpS(useJ),ANNdecisionData.correctedTable.ANNbpS(useJ),'.k',...
                'DisplayName',['corrected decisions (r \approx ',num2str(rhoc,3),')'],'MarkerSize',24)
        else
            [rhoc,pvalc] = corr(ANNdecisionData.correctedTable.EUCASTbpS(useJ),ANNdecisionData.correctedTable.ANNbpS(useJ),'type','Pearson');
            plot(ANNdecisionData.correctedTable.EUCASTbpS(useJ),ANNdecisionData.correctedTable.ANNbpS(useJ),'.k',...
                'DisplayName',['corrected decisions (r \approx ',num2str(rhoc,3),')'],'MarkerSize',24)
        end
    end
    X = -5:0.1:7;
    plot(X,X,'-k','DisplayName','x=y')
    plot(X,mdlS.feval(X),'--k','DisplayName','linear regression')
    plot(X,NLmdlS.feval(X),':k','DisplayName','cubic regression')
    %set(gca,'Xtick',P.lDosages)
    %set(gca,'XtickLabels',P.lDosages)
    if CEflag
        xlabel('CLSI breakpoint S')
    else
        xlabel('EUCAST breakpoint S')
    end
    ylabel('Neural Network breakpoint S')
    axis tight
    legend('location','northwest')

    figure(figNum+1)
    if CEflag
        [rho,pval] = corr(ANNdecisionData.table.CLSIbpR(useJ),ANNdecisionData.table.ANNbpR(useJ),'type','Pearson');
        NLmdlR = fitnlm(ANNdecisionData.table.CLSIbpR(useJ),ANNdecisionData.table.ANNbpR(useJ),f,[0 1 0 2])
        mdlR = fitlm(ANNdecisionData.table.CLSIbpR(useJ),ANNdecisionData.table.ANNbpR(useJ));
    else
        [rho,pval] = corr(ANNdecisionData.table.EUCASTbpR(useJ),ANNdecisionData.table.ANNbpR(useJ),'type','Pearson');
        NLmdlR = fitnlm(ANNdecisionData.table.EUCASTbpR(useJ),ANNdecisionData.table.ANNbpR(useJ),f,[0 1 0 2])
        mdlR = fitlm(ANNdecisionData.table.EUCASTbpR(useJ),ANNdecisionData.table.ANNbpR(useJ));
    end

    if CEflag
        plot(ANNdecisionData.table.CLSIbpR(useJ),ANNdecisionData.table.ANNbpR(useJ),'o','color',dotColour,'DisplayName',['raw decisions (r \approx ',num2str(rho,3),')'],'MarkerSize',14)
    else
        plot(ANNdecisionData.table.EUCASTbpR(useJ),ANNdecisionData.table.ANNbpR(useJ),'o','color',dotColour,'DisplayName',['raw decisions (r \approx ',num2str(rho,3),')'],'MarkerSize',14)
    end
    hold on
    try
        if CEflag
            [rhoc,pvalc] = corr(ANNdecisionData.correctedTable.CLSIbpR(useJ),ANNdecisionData.correctedTable.ANNbpR(useJ),'type','Pearson');
            plot(ANNdecisionData.correctedTable.CLSIbpR(useJ),ANNdecisionData.correctedTable.ANNbpR(useJ),'.k',...
                'DisplayName',['corrected decisions (r \approx ',num2str(rhoc,3),')'],'MarkerSize',24)
        else
            [rhoc,pvalc] = corr(ANNdecisionData.correctedTable.EUCASTbpR(useJ),ANNdecisionData.correctedTable.ANNbpR(useJ),'type','Pearson');
            plot(ANNdecisionData.correctedTable.EUCASTbpR(useJ),ANNdecisionData.correctedTable.ANNbpR(useJ),'.k',...
                'DisplayName',['corrected decisions (r \approx ',num2str(rhoc,3),')'],'MarkerSize',24)
        end
    end
    X = -5:0.1:10;
    plot(X,X,'-k','DisplayName','x=y')
    plot(X,mdlR.feval(X),'--k','DisplayName','linear regression')
    plot(X,NLmdlR.feval(X),':k','DisplayName','cubic regression')
    if CEflag
        xlabel('CLSI breakpoint R')
    else
        xlabel('EUCAST breakpoint R')
    end
    ylabel('Neural Network breakpoint R')
    axis tight
    legend('location','northwest')

    %set(gca,'Xtick',P.lDosages)
    %set(gca,'XtickLabels',P.lDosages)

end