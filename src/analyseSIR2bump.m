function [ignoreListOut,mdl,tbl,targetString] = analyseSIR2bump(SIRData,J,ignoreList)
    
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
    GMMRmean = data.GMMRmean;
    GMMRsigma = data.GMMRsigma;
    SRboundary = data.SRboundary;
    CLSIsticS = data.CLSIsticS;
    CLSIsticR = data.CLSIsticR;
    eucastBPS = data.eucastBPS;
    eucastBPR = data.eucastBPR;
    Tecoff = data.Tecoff;
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

    for j = 1:8
        figure(j)
        set(j,'pos',[90         882        600         455])
    end

    %The 2-bump stuff:
    F = (Nclusters == 2);
    F(ignoreList) = 0;

    F = F & GMMconverged;
    %F = F & ~isnan(SRboundary);
    F = F & ~isnan(targetVector);
    
    Fex = F & ~isnan(Xecoff);
    Fet = F & ~isnan(Tecoff);

    Fetex = Fet & Fex;

    figure(1)
    plot(GMMRsquared(Fex),XecoffCorrelation(Fex),'.k','markersize',16,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName','x=y')
    xlim(xL)
    xlabel('2-normal (GMM) fit R^2')
    ylabel('est. ECOFF correlation')
    legend('location','northwest')

    figure(2)
    rho = corr(Xecoff(Fex),SRboundary(Fex));
    plot(Xecoff(Fex),SRboundary(Fex),'.k','markersize',14,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('SR boundary')
    xlabel('est. ECOFF')
    legend('location','northwest');
    %,'HandleVisibility','off'

    figure(22)
    rho = corr(Tecoff(Fex),SRboundary(Fex));
    plot(Tecoff(Fex),SRboundary(Fex),'.k','markersize',14,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('SR boundary')
    xlabel('(T)ECOFF')
    legend('location','northwest');
    %,'HandleVisibility','off'

    figure(3)
    rho = corr(Xecoff(Fex),GMMSmean(Fex));
    plot(Xecoff(Fex),GMMSmean(Fex),'.k','markersize',14,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('S mean')
    xlabel('est. ECOFF')
    legend('location','northwest');
    %,'HandleVisibility','off'

    figure(33)
    rho = corr(Tecoff(Fex),GMMSmean(Fex) + 2*GMMSsigma(Fex));
    plot(Tecoff(Fex),GMMSmean(Fex) + 2*GMMSsigma(Fex),'.k','markersize',14,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('S mean + 2\sigma')
    xlabel('(T)ECOFF')
    legend('location','northwest');
    %,'HandleVisibility','off'
    
    figure(4)
    rho = corr(Xecoff(Fex),GMMRmean(Fex));
    plot(Xecoff(Fex),GMMRmean(Fex),'.k','markersize',14,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('R mean')
    xlabel('est. ECOFF')
    legend('location','northwest');
    %,'HandleVisibility','off'

    figure(44)
    rho = corr(Tecoff(Fex),GMMRmean(Fex));
    plot(Tecoff(Fex),GMMRmean(Fex),'.k','markersize',14,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('R mean')
    xlabel('(T)ECOFF')
    legend('location','northwest');
    %,'HandleVisibility','off'
    
    figure(5)
    rho = corr(Xecoff(Fex),SRboundary(Fex));
    plot(Xecoff(Fex),SRboundary(Fex),'.k','markersize',14,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('SR boundary')
    xlabel('est. ECOFF')
    legend('location','northwest');
    %,'HandleVisibility','off'

    figure(55)
    rho = corr(Tecoff(Fex),SRboundary(Fex));
    plot(Tecoff(Fex),SRboundary(Fex),'.k','markersize',14,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('SR boundary')
    xlabel('(T)ECOFF')
    legend('location','northwest');
    %,'HandleVisibility','off'
    
    figure(6)
    rho = corr(targetVector(Fet),Tecoff(Fet));    
    plot(targetVector(Fet),Tecoff(Fet),'.k','markersize',14,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fet)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('(T)ECOFF')
    xlabel(targetString)
    legend('location','northwest');
    %,'HandleVisibility','off'

    tbl = table(GMMSmean(F),sqrt(GMMSsigma(F)),GMMRmean(F),sqrt(GMMRsigma(F)),SRboundary(F),targetVector(F),...
        'VariableNames',{'Smean','Sstd','Rmean','Rstd','SRboundary',targetStringMDL});

    %use variance instead of S.D.:
    %tbl = table(GMMSmean(F),(GMMSsigma(F)),GMMRmean(F),(GMMRsigma(F)),SRboundary(F),targetVector(F),...
    %    'VariableNames',{'Smean','Sstd','Rmean','Rstd','SRboundary',targetStringMDL});

    mdl = stepwiselm(tbl,[targetStringMDL,'~1+Smean+Rmean+Sstd+Rstd+SRboundary'],...
        'ResponseVar',targetStringMDL,'PredictorVars',{'Smean','Sstd','Rmean','Rstd','SRboundary'},...
        'Upper','poly11111');
    disp(mdl)

    figure(7)
    Y = mdl.Fitted;
    tVF = targetVector(F);
    rho = corr(tVF,Y);    
    plot(tVF,Y,'.k','markersize',28,'DisplayName',['2-cluster datasets (N=',num2str(sum(F)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    G = (abs(Y-tVF) > pid);
    plot(tVF(G),Y(G),'or','markersize',12,'DisplayName',['>',num2str(pid),' MIC difference (',num2str(sum(G)),' datapoints)'])

    ylabel('linear model prediction')
    xlabel(targetString)
    title('prediction from stepwise linear model')
    legend('location','northwest');

    figure(8)
    rho = corr(Tecoff(Fetex),Xecoff(Fetex));    
    plot(Tecoff(Fetex),Xecoff(Fetex),'.k','markersize',14,'DisplayName',['2-cluster datasets (N=',num2str(sum(Fetex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    plot([xL(1) xL(2)],[pid+xL(1) pid+xL(2)],'--k','linewidth',1,'DisplayName','discrepancy limits')
    plot([xL(1) xL(2)],[-pid+xL(1) -pid+xL(2)],'--k','linewidth',1,'HandleVisibility','off')
    ylabel('est. ECOFF')
    xlabel('(T)ECOFF')
    legend('location','northwest');

    G = (abs(Tecoff-Xecoff)>pid);
    G = G & Fetex;
    ignoreListOut = find(G);

    plot(Tecoff(G),Xecoff(G),'or','markersize',12,'DisplayName',['>',num2str(pid),'MIC difference'])
    if not(isempty(ignoreListOut))
        disp('-----------------------')
        disp('Discrepancies found ...')
        disp('-----------------------')
    end
    for i = 1:length(ignoreListOut)
        u = Tecoff(ignoreListOut(i));
        v = Xecoff(ignoreListOut(i));
        text(u,v,{drugs{ignoreListOut(i)},bugs{ignoreListOut(i)}});
        disp([drugs{ignoreListOut(i)},' & ',bugs{ignoreListOut(i)},' at input SIRData file address ',num2str(2+ignoreListOut(i))])
    end
    if not(isempty(ignoreListOut))
        disp('-----------------------')
    end
    
    if J == 1
        figure(9)
        N = 10;
        Constant = mdl.Coefficients.Estimate(1);
        alpha = mdl.Coefficients.Estimate(2);
        beta = mdl.Coefficients.Estimate(3);
        gamma = mdl.Coefficients.Estimate(3);
        sm = GMMSmean(F);
        smRange = linspace(min(sm),max(sm),N);
        srb = SRboundary(F);
        srbRange = linspace(min(srb),max(srb),N);
        rm = GMMRmean(F);
        rmRange = linspace(min(rm),max(rm),N);
    
        %assume beta = 0
        mdlFunction = @(x,z)(Constant + alpha*x + gamma*z);
        
        [X,Y] = meshgrid(smRange,srbRange);
        surf(X,Y,mdlFunction(X,Y))
        hold on
        plot3(sm,srb,targetVector(F),'.k','MarkerSize',12)
        xlabel('mean S')
        ylabel('SR boundary')
        zlabel(targetString)
        axis tight
        colorbar
    end
end
