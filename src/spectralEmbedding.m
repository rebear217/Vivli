function spectralEmbedding(allEUCASTHistograms,consolidatedSIRData,CEflag,weight)

    if nargin < 3
        CEflag = 1;
    end

    if nargin < 4
        weight = ones(size((lDosages)));
    end
    weight = weight/sum(weight);

    close all

    data = extractSIRdataColumns(consolidatedSIRData);
    lDosages = data.lDosages;

    %lDosagesString = data.lDosagesString;
    bpSclsi = data.CLSIsticS;
    bpRclsi = data.CLSIsticR;
    bpSeucast = data.eucastBPS;
    bpReucast = data.eucastBPR;
    drugs = data.drugs;
    bugs = data.bugs;
    Nclusters = 25;

    clear data

    if CEflag
        breakpointsS = bpSclsi;
        breakpointsR = bpRclsi;
    else
        breakpointsS = bpSeucast;
        breakpointsR = bpReucast;
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
    N = length(usePositions);

    Wmatrix = zeros(N,N);
    bpDiffR = zeros(N,N);
    bpDiffS = zeros(N,N);
    for j = 1:N
        J = histoPositions(usePositions(j));
        HJ = allHistograms(J,:);
        bpJS = breakpointsS(usePositions(j));
        bpJR = breakpointsR(usePositions(j));
        for i = (j+1):N
            I = histoPositions(usePositions(i));
            bpIS = breakpointsS(usePositions(i));
            bpIR = breakpointsR(usePositions(i));
            HI = allHistograms(I,:);
            Wmatrix(i,j) = myNorm(HI-HJ)*sign(max(HI-HJ));
            Wmatrix(j,i) = myNorm(HI-HJ)*sign(max(HJ-HI));
            bpDiffS(i,j) = bpJS-bpIS;
            bpDiffS(j,i) = bpJS-bpIS;
            bpDiffR(i,j) = bpJR-bpIR;
            bpDiffR(j,i) = bpJR-bpIR;
        end
    end

    %LipMatrixS = bpDiffS ./ Wmatrix;
    %LipMatrixR = bpDiffR ./ Wmatrix;

    %cgW = clustergram(Wmatrix-0.5,'Standardize','none','Symmetric','true','Colormap','hot');
    %%clustergram(bpDiffS,'Standardize','none','Symmetric','true','Colormap','hot');
    %%clustergram(bpDiffR,'Standardize','none','Symmetric','true','Colormap','hot');

    function idx = weightOrdering(allHs)
        specvec = spectralcluster(allHs,Nclusters);
        specSizes = zeros(size(specvec));
        for jf = 1:Nclusters
            theseIndicies = (specvec == jf);
            mvec = mean(allHs(theseIndicies,:),1);
            specSizes(theseIndicies) = sum(lDosages.*mvec)/sum(mvec);
        end
        [~,idx] = sort(specSizes);
    end

    %This usess the method in https://graphics.stanford.edu/courses/cs205a-13-fall/assets/notes/chapter5.pdf
    %to order histograms according to optimally sorting their pairwise differences:
    %it does not work very well:
    %Wvec = sum(Wmatrix,1);
    %A = diag(Wvec);
    %[V,D] = eigs(2*(A-Wmatrix),2,1e-5);
    %vec = V(:,2);

    %[V,D] = eig(2*(A-Wmatrix));
    %d = diag(D);
    %d = d(:);
    %[m,jj] = mink(d,2);
    %vec = V(:,jj(2));

    %vec = vec / norm(vec,1);
    %[~,Iv] = sort(vec);

    figure(1)
    set(1,'pos',[109    84   844   613])
    imagesc(Wmatrix);
    title('EUCAST histograms with breakpoints pairwise norm differences')
    colorbar

    figure(2)
    set(2,'pos',[79         469        1226         548])
    idx = weightOrdering(allHistograms(histoPositions(usePositions),:));
    imagesc(flipud(allHistograms(histoPositions(usePositions(idx)),:)'));
    hold on
    
    for i = 1:N
        p1 = plot([i i],[10-bpSclsi(usePositions(idx(i))) 10-bpRclsi(usePositions(idx(i)))],'-r','linewidth',1);
        if bpSeucast(usePositions(idx(i))) > -7
            p2 = plot([i i],[10-bpSeucast(usePositions(idx(i))) 10-bpReucast(usePositions(idx(i)))],'-b','linewidth',2);
        end
    end
    
    CLSIuse = bpSclsi(usePositions(idx));
    nCLSIuse = sum(~isnan(CLSIuse));
    EUuse = bpSeucast(usePositions(idx));
    nEUuse = sum(~isnan(EUuse));

    p3 = plot(20,-20,'-k');
    
    title('Heatmap of observed frequencies')
    legend([p1,p2,p3],{['CLSI SR breakpoint range (# data ',num2str(nCLSIuse),')'],...
        ['EUCAST SR breakpoint range  (# data ',num2str(nEUuse),')'],'MIC distributions viewed from above'},'location','southeast')
    colormap(1-bone)
    ylabel('log2 MIC dose')
    set(gca,'Ytick',1:19)
    set(gca,'Yticklabel',fliplr(round(lDosages,2)))
    xlabel('Spectrally ordered EUCAST histograms with breakpoints')
    colorbar
    drawnow

    %pause

    W = triu(Wmatrix,1);
    W = W(:);
    newW = W(W>0);
    bdifR = triu(bpDiffR,1);
    bdifR = bdifR(:);
    bdifR = bdifR(W > 0);
    bdifS = triu(bpDiffS,1);
    bdifS = bdifS(:);
    bdifS = bdifS(W > 0);
    W = newW;

    fairnessS = bdifS./W;
    fairnessR = bdifR./W;
    zfS = sum((bdifS == 0));
    zfR = sum((bdifR == 0));
    
    figure(3)
    set(3,'pos',[911   535   699   386])
    histogram(fairnessS)
    %histogram(fairnessS(~(fairnessS == 0)),'DisplayName',[num2str(zfS),' have zero S-sensitivity coefficient (not shown)'])
    xlabel('Sensitivity coefficient (L)')
    text(200,4000,'$L = \frac{b(H)-b(K)}{\|H-K\|_{w}}$','interpreter','latex','FontSize',20)

    hold on
    histogram(fairnessR)
    %histogram(fairnessR(~(fairnessR == 0)),'DisplayName',[num2str(zfR),' have zero R-sensitivity coefficient (not shown)'])
    xlabel('Sensitivity coefficient (L)')
    ylabel('# observations')
    axis tight
    xlim([-500 500])
    legend()

    figure(4)
    set(4,'pos',[222   310   845   356])
    subplot(1,2,1)
    plot(W,bdifS+0.05*randn(size(W)),'.','color',[1 1 1]*0.7);
    ylabel('S breakpoint difference (+noise)')
    xlabel('histogram difference')
    ylim([-10 10])
    hold on
    subplot(1,2,2)
    plot(W,bdifR+0.05*randn(size(W)),'.','color',[1 1 1]*0.7);
    ylabel('R breakpoint difference (+noise)')
    xlabel('histogram difference')
    ylim([-10 10])
    hold on

    %zfSI = find(bdifS == 0);
    %zfRI = find(bdifR == 0);

    figure(5)
    set(5,'pos',[350   312   319   677])
    subplot(3,1,1)
    histogram(W(fairnessS == 0),'DisplayName',['equal S bps (',num2str(zfS),' datapoints)'])
    %hold on
    xlabel('histogram difference')
    ylabel('# observations')
    legend()
    ylim([0 1000])

    subplot(3,1,2)
    histogram(W(fairnessR == 0),'DisplayName',['equal R bps (',num2str(zfR),' datapoints)'])
    xlabel('histogram difference')
    ylabel('# observations')
    legend()
    ylim([0 1000])

    zfSR = (fairnessS == 0) & (fairnessR == 0);

    subplot(3,1,3)
    h = histogram(W(zfSR),'DisplayName',['equal S&R bps (',num2str(sum(zfSR)),' datapoints)']);
    xlabel('histogram difference')
    ylabel('# observations')
    legend()
    ylim([0 1000])

    figure(6)
    set(6,'pos',[138     1   544   696])
    %idx = weightOrdering(allHistograms(histoPositions,:));
    %imagesc(allHistograms(histoPositions(idx),:));
    idx = weightOrdering(allHistograms);
    imagesc(allHistograms(idx,:));
    title('All spectrally ordered histograms')
    colormap(1-bone)
    xlabel('log2 MIC dose')
    set(gca,'Xtick',1:19)
    set(gca,'Xticklabel',round(lDosages,2))
    ylabel('heatmap of observed frequencies')
    colorbar
    
    %figure(100)
    %view(cgW);

    panels = 6;
    wzi = find(cumsum(fliplr(h.Values)) > panels^2,1);
    threshold = h.BinEdges(end-wzi);
    [WFzSRi,WFzSRj] = find((Wmatrix > threshold) & (bpDiffS == 0) & (bpDiffR == 0));
    
    figure(4)
    for ii = 1:min(panels^2,length(WFzSRi))
        i = WFzSRi(ii);
        j = WFzSRj(ii);

        subplot(1,2,1)
        plot(Wmatrix(i,j),0,'.r','markersize',14)
        subplot(1,2,2)
        plot(Wmatrix(i,j),0,'.r','markersize',14)
    end

    %figure(7)
    %set(7,'pos',[90     1   991   696])
    for ii = 1:min(panels^2,length(WFzSRi))
        %subplot(panels,panels,ii);
        figure(7+ii)
        i = WFzSRi(ii);
        j = WFzSRj(ii);
        I = histoPositions(usePositions(i));
        HI = allHistograms(I,:);
        J = histoPositions(usePositions(j));
        
        HJ = allHistograms(J,:);
        plot(lDosages,HI,'.-k','linewidth',1,'DisplayName',[drugs{usePositions(i)},' & ',bugs{usePositions(i)}],'markersize',28)
        hold on
        plot(lDosages,HJ,'.-b','linewidth',1,'DisplayName',[drugs{usePositions(j)},' & ',bugs{usePositions(j)}],'markersize',28)
        plot([breakpointsS(usePositions(i)) breakpointsS(usePositions(i))],[0 1],'--k','linewidth',1,'DisplayName','breakpoint S');
        plot([breakpointsR(usePositions(i)) breakpointsR(usePositions(i))],[0 1],'-k','linewidth',1,'DisplayName','breakpoint R');
        %plot([breakpointsS(usePositions(j)) breakpointsS(usePositions(j))],[0 1],'--r','linewidth',1);
        %plot([breakpointsR(usePositions(j)) breakpointsR(usePositions(j))],[0 1],'-r','linewidth',1);
        xlabel('log2 MIC')
        ylabel('frequency')
        axis tight
        ylim([0 1])
        legend();
        text(-8.5,0.1,['$\|H_i-H_j\|_{w} = $',num2str(Wmatrix(i,j),5)],'interpreter','latex')
        title('Same breakpoints, maximally distant histograms')
    end

    wzi = find(cumsum(h.Values) > panels^2,1);
    if wzi <= 1
        wzi = 2;
    end
    threshold = h.BinEdges(wzi);
    [WFzSRi,WFzSRj] = find((Wmatrix < threshold) & (abs(bpDiffS) > 0.5) & (abs(bpDiffR) > 0.5));


    ff = gcf;
    %figure(ff.Number)
    %set(ff.Number,'pos',[90     1   991   696])
    for ii = 1:min(panels^2,length(WFzSRi))
        %subplot(panels,panels,ii);
        figure(ff.Number+ii)
        i = WFzSRi(ii);
        j = WFzSRj(ii);
        I = histoPositions(usePositions(i));
        HI = allHistograms(I,:);
        J = histoPositions(usePositions(j));
        
        HJ = allHistograms(J,:);
        plot(lDosages,HI,'.-k','linewidth',1,'DisplayName',[drugs{usePositions(i)},' & ',bugs{usePositions(i)}],'markersize',28)
        hold on
        plot(lDosages,HJ,'.-b','linewidth',1,'DisplayName',[drugs{usePositions(j)},' & ',bugs{usePositions(j)}],'markersize',28)
        plot([breakpointsS(usePositions(i)) breakpointsS(usePositions(i))],[0 1],'--k','linewidth',1,'DisplayName','1st breakpoint S');
        plot([breakpointsR(usePositions(i)) breakpointsR(usePositions(i))],[0 1],'-k','linewidth',1,'DisplayName','1st breakpoint R');
        plot([breakpointsS(usePositions(j)) breakpointsS(usePositions(j))],[0 1],'--b','linewidth',1,'DisplayName','2nd breakpoint S');
        plot([breakpointsR(usePositions(j)) breakpointsR(usePositions(j))],[0 1],'-b','linewidth',1,'DisplayName','2nd breakpoint R');
        xlabel('log2 MIC')
        ylabel('frequency')
        axis tight
        ylim([0 1])
        set(gca,'XTick',round(lDosages,2))
        set(gca,'XTickLabel',round(lDosages,2))
        legend();
        text(-8.5,0.1,['$\|H_i-H_j\|_{w} = $',num2str(Wmatrix(i,j),5)],'interpreter','latex')
        title('Different breakpoints, nearby histograms')
    end

    figure(4)
    for ii = 1:min(panels^2,length(WFzSRi))
        i = WFzSRi(ii);
        j = WFzSRj(ii);

        subplot(1,2,1)
        plot(Wmatrix(i,j),bpDiffS(i,j),'.b','markersize',18)
        subplot(1,2,2)
        plot(Wmatrix(i,j),bpDiffR(i,j),'.b','markersize',18)
    end
    
    subplot(1,2,1)
    xL1 = xlim;
    subplot(1,2,2)
    xL2 = xlim;

    %{
    Q = -0.1:0.01:0.1;
    y = 2:0.01:10;
    for j = 1:length(Q)
        q = Q(j);
        xp = q*(y + sqrt(y.*2 - 4))/2;
        xm = q*(y - sqrt(y.*2 - 4))/2;

        s = (j-1)/(length(Q)-1);
        colr = [1-s 0 s];

        subplot(1,2,1)
        plot(xp,y,'-k','Linewidth',1,'color',colr);
        plot(xm,y,'-k','Linewidth',1,'color',colr);
        plot(y,xp,'-k','Linewidth',1,'color',colr);
        plot(y,xm,'-k','Linewidth',1,'color',colr);
        subplot(1,2,2)
        plot(-xp,y,'-k','Linewidth',1,'color',colr);
        plot(-xm,y,'-k','Linewidth',1,'color',colr);
        plot(y,xp,'-k','Linewidth',1,'color',colr);
        plot(y,xm,'-k','Linewidth',1,'color',colr);
    end
    %}
    subplot(1,2,1)
    %ezcontour('-2+x./y + y./x',[0,0.1],[-10 10],100)
    %ezcontour('y./x',[0,0.1],[-10 10],100)
    xlim([0 xL1(2)]);
    subplot(1,2,2)
    xlim([0 xL2(2)]);

    function y=myNorm(v)
        %N = length(v);
        %w = (1:N);
        %w = w/N;
        y = norm(weight.*v);
        %y = norm(v);
    end

end