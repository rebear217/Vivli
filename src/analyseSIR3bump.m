function [ignoreListOut,mdl,tbl,targetString] = analyseSIR3bump(SIRData,J,ignoreList)

    Params = globalParameterValues();
    pid = Params.ignoreDiscrepancy;

    if nargin < 2
        J = 1;
    end
    if nargin < 3
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
    eucastBPS = data.eucastBPS;
    eucastBPR = data.eucastBPR;
    Xecoff = data.Xecoff;
    XecoffCorrelation = data.XecoffCorrelation;

    clear data
    switch J
        case 1
            targetVector = Tecoff;
            targetString = '(T)ECOFF';
            targetStringMDL = 'TECOFF';
        case 2
            targetVector = CLSIsticS;
            targetString = 'CLSISTICS';
            targetStringMDL = 'CLSISTICS';
        case 3
            targetVector = CLSIsticR;
            targetString = 'CLSISTICR';
            targetStringMDL = 'CLSISTICR';
        case 4
            targetVector = eucastBPS;
            targetString = 'eucastBPS';
            targetStringMDL = 'eucastBPS';
        case 5
            targetVector = eucastBPR;
            targetString = 'eucastBPR';
            targetStringMDL = 'eucastBPR';
    end

    for j = 1:4
        figure(j)
        set(j,'pos',[90         882        600         455])
    end

    %The 3-bump stuff:
    F = (Nclusters == 3);
    F(ignoreList) = 0;
    F = F & GMMconverged;
    %F = F & ~isnan(SRboundary);
    %F = F & ~isnan(SIboundary);
    %F = F & ~isnan(IRboundary);
    F = F & ~isnan(targetVector);    
    
    Fex = F & ~isnan(Xecoff);
    Fet = F & ~isnan(Tecoff);

    Fetex = Fet & Fex;

    figure(1)
    rho = corr(Tecoff(Fet),targetVector(Fet));    
    plot(Tecoff(Fet),targetVector(Fet),'.k','markersize',14,'DisplayName',['3-cluster datasets (N=',num2str(sum(Fet)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    xlabel('TECOFF')
    ylabel(targetString)
    legend('location','northwest');

    figure(2)
    rho = corr(Xecoff(Fex),targetVector(Fex));    
    plot(Xecoff(Fex),targetVector(Fex),'.k','markersize',14,'DisplayName',['3-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    xlabel('est. ECOFF')
    ylabel(targetString)
    legend('location','northwest');

    figure(3)
    rho = corr(Xecoff(Fetex),Tecoff(Fetex));    
    plot(Xecoff(Fetex),Tecoff(Fetex),'.k','markersize',14,'DisplayName',['3-cluster datasets (N=',num2str(sum(Fetex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    plot([xL(1) xL(2)],[pid+xL(1) pid+xL(2)],'--k','linewidth',1,'DisplayName','discrepancy limits')
    plot([xL(1) xL(2)],[-pid+xL(1) -pid+xL(2)],'--k','linewidth',1,'HandleVisibility','off')

    ylabel('est. ECOFF')
    xlabel('TECOFF')
    legend('location','northwest');    

    G = (abs(Tecoff-Xecoff)>pid);
    G = G & Fetex;
    ignoreListOut = find(G);

    plot(Xecoff(G),Tecoff(G),'or','markersize',12,'DisplayName',['>',num2str(pid),'MIC difference'])
    if not(isempty(ignoreListOut))
        disp('-----------------------')
        disp('Figure 3: Large MIC differences found for ...')
        disp('-----------------------')
    end
    for i = 1:length(ignoreListOut)
        u = Tecoff(ignoreListOut(i));
        v = Xecoff(ignoreListOut(i));
        text(v,u,{drugs{ignoreListOut(i)},bugs{ignoreListOut(i)}});
        disp([drugs{ignoreListOut(i)},' & ',bugs{ignoreListOut(i)},' at input SIRData file address ',num2str(2+ignoreListOut(i))])
    end
    if not(isempty(ignoreListOut))
        disp('-----------------------')
    end
    
    tbl = table(GMMSmean(F),sqrt(GMMSsigma(F)),GMMImean(F),sqrt(GMMIsigma(F)),...
                GMMRmean(F),sqrt(GMMRsigma(F)),SIboundary(F),IRboundary(F),SRboundary(F),targetVector(F),...
                'VariableNames',{'Smean','Sstd','Imean','Istd','Rmean','Rstd','SIboundary','IRboundary',...
                'SRboundary',targetStringMDL});

    mdl = stepwiselm(tbl,[targetStringMDL,'~1+Smean+Imean+Rmean+Sstd+Istd+Rstd+SIboundary+IRboundary+SRboundary'],...
        'ResponseVar',targetStringMDL,'PredictorVars',...
        {'Smean','Sstd','Imean','Istd','Rmean','Rstd','SIboundary','IRboundary','SRboundary'},...
        'Upper','poly111111111');
    disp(mdl)

    figure(4)
    tVF = targetVector(F);
    rho = corr(tVF,mdl.Fitted);    
    plot(tVF,mdl.Fitted,'.k','markersize',28,'DisplayName',['3-cluster datasets (N=',num2str(sum(F)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('linear model prediction')
    xlabel(targetString)
    title('prediction from stepwise linear model')
    legend('location','northwest');
    G = (abs(mdl.Fitted-tVF) > pid);
    plot(tVF(G),mdl.Fitted(G),'or','markersize',12,'DisplayName',['>',num2str(pid),' MIC difference (',num2str(sum(G)),' datapoints)'])
    
end