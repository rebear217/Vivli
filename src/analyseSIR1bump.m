function ignoreListOut = analyseSIR1bump(SIRData,ignoreList)

    Params = globalParameterValues();
    pid = Params.ignoreDiscrepancy;

    if nargin < 2
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
    %CLSIsticL = data.CLSIsticL;
    %CLSIsticR = data.CLSIsticR;
    Tecoff = data.Tecoff;
    Xecoff = data.Xecoff;

    clear data

    %where the Tecoff might be unacceptably low??:
    %sum(Tecoff < GMMSmean-GMMSsigma)
    %sum(Tecoff < Xecoff)

    %The 1-bump stuff:
    for j = 1:5
        figure(j)
        set(j,'pos',[90         882        600         455])
    end

    F = (Nclusters == 1);
    F(ignoreList) = 0;

    F = F & GMMconverged;

    Fex = F & ~isnan(Xecoff);
    Fet = F & ~isnan(Tecoff);

    Fetex = Fet & Fex;

    ECpredictor = (Xecoff(Fex) - GMMSmean(Fex))./sqrt(GMMSsigma(Fex));
    muECm = mean(ECpredictor);

    figure(1)
    plot(GMMRsquared(Fex),ECpredictor,'.k','markersize',16,'DisplayName',['1-cluster datasets (N=',num2str(sum(Fex)),')'])
    axis tight
    xL = xlim;
    xL = [xL(1) 1];
    xlim(xL)
    hold on
    plot([xL(1) xL(2)],[muECm muECm],'-k','linewidth',1,'DisplayName','mean ECOFF predictor (p)')
    xlabel('normal (GMM) fit R^2')
    ylabel('ECOFF predictor (p_i)')
    legend('location','northwest')

    figure(2)
    histogram(ECpredictor,100)
    ylabel('number observed')
    xlabel('ECOFF predictor (p_i)')

    figure(3)
    Y = GMMSmean + muECm*sqrt(GMMSsigma);
    rho = corr(Xecoff(Fex),Y(Fex));
    plot(Xecoff(Fex),Y(Fex),'.k','markersize',16,'DisplayName',['1-cluster datasets (N=',num2str(sum(Fex)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    ylabel('normal predicted ECOFF')
    xlabel('est. ECOFF')
    text(0,-4,'ECOFF${}_i = \mu_i + p\sqrt\sigma_i$','interpreter','latex')
    text(0,-5,['$p\approx',num2str(muECm),'$'],'interpreter','latex')
    legend('location','northwest')

    figure(4)
    rho = corr(Tecoff(Fet),Y(Fet));
    plot(Tecoff(Fet),Y(Fet),'.k','markersize',28,'DisplayName',['1-cluster datasets (N=',num2str(sum(Fet)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    %plot([xL(1) xL(2)],[pid+xL(1) pid+xL(2)],'--k','linewidth',1,'DisplayName','discrepancy limits')
    %plot([xL(1) xL(2)],[-pid+xL(1) -pid+xL(2)],'--k','linewidth',1,'HandleVisibility','off')
    ylabel('normal predicted (T)ECOFF')
    xlabel('(T)ECOFF')
    legend('location','northwest');
    %,'HandleVisibility','off'
    G = (abs(Tecoff-Y) > pid);
    G = G & Fet;
    ignoreListOut = find(G);
    text(0,-4,'ECOFF${}_i = \mu_i + p\sqrt\sigma_i$','interpreter','latex')
    text(0,-5,['$p\approx',num2str(muECm),'$'],'interpreter','latex')

    plot(Tecoff(G),Y(G),'or','markersize',12,'DisplayName',['>',num2str(pid),' MIC difference (',num2str(sum(G)),' datapoints)'])
    if not(isempty(ignoreListOut))
        disp('-----------------------')
        disp('Model-data differences found for ...')
        disp('-----------------------')
    end
    for i = 1:length(ignoreListOut)
        u = Tecoff(ignoreListOut(i));
        v = Y(ignoreListOut(i));
        text(u-1,v-0.5,[drugs{ignoreListOut(i)},' & ',bugs{ignoreListOut(i)}]);
        disp([drugs{ignoreListOut(i)},' & ',bugs{ignoreListOut(i)},' at input SIRData file address ',num2str(2+ignoreListOut(i))])
    end
    if not(isempty(ignoreListOut))
        disp('-----------------------')
    end

    figure(5)
    rho = corr(Tecoff(Fetex),Xecoff(Fetex));    
    plot(Tecoff(Fetex),Xecoff(Fetex),'.k','markersize',14,'DisplayName',['1-cluster datasets (N=',num2str(sum(Fetex)),')'])
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
    %,'HandleVisibility','off'

    G = (abs(Tecoff-Xecoff) > pid);
    G = G & Fetex;
    fg = find(G);
    plot(Tecoff(G),Xecoff(G),'or','markersize',12,'DisplayName',['>',num2str(pid),' MIC difference (',num2str(sum(G)),' datapoints)'])
    for i = 1:length(fg)
        u = Tecoff(fg(i));
        v = Xecoff(fg(i));
        text(u-4,v,{drugs{fg(i)},bugs{fg(i)}});
    end

end