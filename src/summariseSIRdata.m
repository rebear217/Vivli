function summariseSIRdata(SIRData,CEflag)

    if nargin < 2
        CEflag = 1;%this is for the use of CLSI bps
        %CEflag = 0 uses EUCAST bps
    end
    
    close all

    ms = 46;
    P = globalParameterValues();
    pid = P.ignoreDiscrepancy;

    red = [1 0 0]*0.7;
    green = [0 1 0]*0.7;
    blue = [0 0 1]*0.7;
    
    figure(1);
    set(1,'pos',[46           1        1235         696]);
    figure(2);
    set(2,'pos',[46           1        1235         696]);
    figure(3);
    set(3,'pos',[46           1        1235         696]);
    figure(4);
    set(4,'pos',[46           1        1235         696]);

    data = extractSIRdataColumns(SIRData);

    if CEflag
        bpcheck = ~isnan(data.CLSIsticR);
    else
        bpcheck = ~isnan(data.eucastBPR);
    end

    TecoffCIhi = data.TecoffCIhi(bpcheck);
    TecoffCIlow = data.TecoffCIlow(bpcheck);

    CLSIsticR = data.CLSIsticR(bpcheck);
    CLSIsticS = data.CLSIsticS(bpcheck);
    eucastBPS = data.eucastBPS(bpcheck);
    eucastBPR = data.eucastBPR(bpcheck);

    Tecoff = data.Tecoff(bpcheck);
    Xecoff = data.Xecoff(bpcheck);
    
    GMMSmean = data.GMMSmean(bpcheck);
    GMMImean = data.GMMImean(bpcheck);
    GMMRmean = data.GMMRmean(bpcheck);
    SIboundary = data.SIboundary(bpcheck);
    SRboundary = data.SRboundary(bpcheck);
    IRboundary = data.IRboundary(bpcheck);

    bugs = data.bugs(bpcheck);
    drugs = data.drugs(bpcheck);
    
    PAnames = cell(length(bugs),1);
    longPAnames = cell(length(bugs),1);
    for j = 1:length(bugs)
        PAnames{j} = shortPAname(drugs{j},bugs{j});
        longPAnames{j} = [drugs{j},'&',bugs{j}];
    end

    %starthere = 50;
    %disp(PAnames{starthere})
    
    corrDatasets = {Tecoff,Xecoff,GMMSmean,GMMImean,GMMRmean,SIboundary,SRboundary,IRboundary};
    corrNames = {'TECOFF','estECOFF','S','I','R','SI','SR','IR'};
    ncd = length(corrDatasets);
    c = eye(ncd);
    for I = 1:ncd
        X = corrDatasets{I};
        NotIsNaNX = ~isnan(X);
        for J = (I+1):ncd
            Y = corrDatasets{J};
            IJ = NotIsNaNX & ~isnan(Y);
            c(I,J) = corr(X(IJ),Y(IJ),'type','Spearman');
            c(J,I) = c(I,J);
        end
    end

    disp(c)
    figure(6)
    imagesc(c);
    shift = zeros(1,ncd)-0.1;
    shift(1) = -0.4;
    shift(2) = -0.4;
    
    for I = 1:ncd
        text(I+shift(I),I,corrNames{I});
    end
    colorbar;
    title('pairwise correlations')
    
    N = length(Tecoff);
    records = 1:N;
    
    TCrange = TecoffCIhi - TecoffCIlow;
    if CEflag
        bpR = CLSIsticR;
        bpS = CLSIsticS;
        BPrange = (CLSIsticR - CLSIsticS)/2;
        BPmid = (CLSIsticR + CLSIsticS)/2;
    else
        bpR = eucastBPR;
        bpS = eucastBPS;
        BPrange = (eucastBPR - eucastBPS)/2;
        BPmid = (eucastBPR + eucastBPS)/2;
    end

    %[CIsorted,CI] = sort(TCrange);
    [CIsorted,CI] = sort(BPmid);
    %[CIsorted,CI] = sort(GMMRmean);    
    
    myXLIM = [149 200];

    for j = 1:2
        figure(j)
        errorbar(records,BPmid(CI(records)),BPrange(CI(records)),'ok','linewidth',1,'markersize',14,'DisplayName','CLSI SR range')
        hold on
        plot(records,GMMSmean(CI(records)),'.','color',green,'markersize',ms,'DisplayName','S mean')
        plot(records,GMMImean(CI(records)),'o','color',blue,'markersize',round(ms/4),'DisplayName','I mean')
        plot(records,GMMRmean(CI(records)),'.','color',red,'markersize',ms,'DisplayName','R mean')
        legend()
        %set(gca,'Xticklabel',data.bugs(CI(records)))
        set(gca,'Ytick',round(P.lDosages,2))
        if j == 2
            set(gca,'Xticklabel',longPAnames(CI(records)))            
            set(gca,'Xtick',1:N)
            xlim(myXLIM)
            xtickangle(45)
        else
            axis tight
        end
        ylim([-9 10])
        ylabel('MIC log2 ug/mL')
    end

    figure(3)
    errorbar(records,Tecoff(CI(records)),TCrange(CI(records)),'ok','linewidth',1,'markersize',14,'DisplayName','TECOFF range')
    hold on
    errorbar(records,BPmid(CI(records)),BPrange(CI(records)),'or','linewidth',1,'markersize',14,'DisplayName','CLSI SR range')
    plot(records,Xecoff(CI(records)),'.','color',[1 1 1]/2,'markersize',20,'DisplayName','est.ECOFF')
    legend()
    %set(gca,'Xticklabel',data.bugs(CI(records)))
    set(gca,'Xticklabel',PAnames(CI(records)))
    set(gca,'Xtick',1:N)
    %xlim([1 100])
    set(gca,'Ytick',round(P.lDosages,2))
    ylim([-10 10])
    ylabel('log2 (ug/mL)')

    figure(4)
    errorbar(records,BPmid(CI(records)),BPrange(CI(records)),'ok','linewidth',1,'markersize',14,'DisplayName','CLSI SR range')
    hold on
    plot(records,SIboundary(CI(records)),'.','color',green,'markersize',ms,'DisplayName','SI boundary')
    plot(records,SRboundary(CI(records)),'o','color',blue,'markersize',round(ms/4),'DisplayName','SR boundary')
    plot(records,IRboundary(CI(records)),'.','color',red,'markersize',ms,'DisplayName','IR boundary')
    legend()
    %set(gca,'Xticklabel',bugs(CI(records)))
    set(gca,'Xticklabel',longPAnames(CI(records)))
    set(gca,'Xtick',1:N)
    %xlim([1 100])
    set(gca,'Ytick',round(P.lDosages,2))
    ylim([-10 10])
    xlim(myXLIM)
    ylabel('log2 (ug/mL)')
    xtickangle(45)

    figure(5)
    Tecoff = data.Tecoff;
    Xecoff = data.Xecoff;
    inT = ~isnan(Tecoff);
    inX = ~isnan(Xecoff);
    T = Tecoff(inT&inX);
    X = Xecoff(inT&inX);
    rho = corr(T,X);
    plot(T,X,'.k','linewidth',1,'markersize',28,'DisplayName',['ECOFF data (N = ',num2str(length(X)),')'])
    hold on
    axis tight
    xL = xlim;
    yL = ylim;
    plot([xL(1) xL(2)],[xL(1) xL(2)],'-k','linewidth',2,'DisplayName',['x=y (\rho\approx',num2str(rho),')'])
    plot([xL(1) xL(2)],[pid+xL(1) pid+xL(2)],'--k','linewidth',1,'HandleVisibility','off')
    plot([xL(1) xL(2)],[-pid+xL(1) -pid+xL(2)],'--k','linewidth',1,'HandleVisibility','off')
    
    G = (abs(T-X) > pid);
    plot(T(G),X(G),'or','markersize',12,'DisplayName',['>',num2str(pid),' MIC difference (',num2str(sum(G)),' datapoints)'])

    %plot([-10 10],[-10 10],'-k','linewidth',1,'DisplayName','equality');
    %text(-9,9,['$\rho\approx$ ',num2str(rho,3)],'interpreter','latex');
    xlabel('published (T)ECOFF')
    ylabel('algorithmic ECOFF')
    legend('location','southeast')

    
    figure(6);
    DataS = GMMSmean - BPmid;
    DataI = GMMImean - BPmid;
    DataR = GMMRmean - BPmid;
    DataBPR = bpR - BPmid;
    DataBPS = bpS - BPmid;

    bins = -9:0.5:9;
    histogram(DataBPR,bins,'FaceColor','k','DisplayName','normalised R breakpoints')
    hold on
    histogram(DataBPS,bins,'FaceColor',[1 1 1]/2,'DisplayName','normalised S breakpoints')
    histogram(DataS,bins,'EdgeColor','g','DisplayName','S cluster centres','DisplayStyle','stairs','linewidth',2)
    histogram(DataR,bins,'EdgeColor','r','DisplayName','R cluster centres','DisplayStyle','stairs','linewidth',2)
    histogram(DataI,bins,'EdgeColor','b','DisplayName','I cluster centres','DisplayStyle','stairs','linewidth',2)
    set(gca,'XtickLabels',-9:9)
    set(gca,'Xtick',-9:9)
    xlabel('log2 \mug/mL')
    ylabel('#observations')
    legend()

    %figure(7)
    %GMMSmean(isnan(GMMSmean)) = 0;
    %GMMImean(isnan(GMMImean)) = 0;
    %GMMRmean(isnan(GMMRmean)) = 0;
    %plot3(GMMSmean,GMMImean,GMMRmean,'ok')

end