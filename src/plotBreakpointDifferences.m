function dtable = plotBreakpointDifferences(bothBreakpointSets)

    close all

    P = globalParameterValues();
    lDosages = P.lDosages;

    figure(1)
    XS = bothBreakpointSets.S_clsi;
    YS = bothBreakpointSets.S_eucast;
    iS9 = find(YS > -9);
    %XS = XS(iS9);
    %YS = YS(iS9);
    
    
    XR = bothBreakpointSets.R_clsi;
    YR = bothBreakpointSets.R_eucast;
    iR9 = find(YR > -9);
    %XR = XR(iR9);
    %YR = YR(iR9);

    grid on

    diffS = YS - XS;
    diffR = YR - XR;
    totDiff = abs(diffS) + abs(diffR);
    [~,iTD] = sort(totDiff,'descend');
    
    plot(diffS,diffR,'.k','markersize',24)
    xlabel('EUCAST - CLSI (S breakpoints)')
    ylabel('EUCAST - CLSI (R breakpoints)')

    fi = find(abs(diffS) > 5 | abs(diffR) > 5);
    hold on
    plot(diffS(fi),diffR(fi),'or','markersize',10)
    for j = 1:length(fi)
        text(diffS(fi(j)),diffR(fi(j)),bothBreakpointSets.FullName_clsi(fi(j)))
    end

    count = 1;
    for j = 1:length(iTD)
        if ~isnan(diffS(iTD(j))) && ~isnan(diffR(iTD(j))) && ismember(iTD(j),iS9) && ismember(iTD(j),iR9)
            bpSclsi(count) = XS(iTD(j));
            bpRclsi(count) = XR(iTD(j));
            bpSeucast(count) = YS(iTD(j));
            bpReucast(count) = YR(iTD(j));
            ds(count) = diffS(iTD(j));
            dr(count) = diffR(iTD(j));
            bugs(count) = bothBreakpointSets.Bacteria_clsi(iTD(j));
            drugs(count) = bothBreakpointSets.Antibiotic_clsi(iTD(j));
            count = count + 1;
        end
    end
    
    dtable = table(ds(:),dr(:),bugs(:),drugs(:),bpSclsi(:),bpRclsi(:),bpSeucast(:),bpReucast(:),...
        'VariableNames',{'dS','dR','bugs','drugs','bpSclsi','bpRclsi','bpSeucast','bpReucast'});
    
    figure(2);
    set(2,'pos',[426         215        1573         450]);

    Drugs = unique(drugs);
    mds = zeros(1,length(Drugs));
    sds = zeros(1,length(Drugs));
    mdr = zeros(1,length(Drugs));
    sdr = zeros(1,length(Drugs));
    
    for j = 1:length(Drugs)
        fj = find(Drugs(j) == drugs);
        mds(j) = mean(ds(fj));
        mdr(j) = mean(dr(fj));
        sds(j) = ste(ds(fj));
        sdr(j) = ste(dr(fj));
    end

    [~,srtmdsj] = sort(mds);
    [~,srtmdrj] = sort(mdr);

    for J = 1:length(Drugs)
        j = srtmdrj(J);

        if J == 1
            errorbar(J-1/4,mds(j),2*sds(j),'.b','markersize',28,'DisplayName','S breakpoint difference mean \pm 2x s.e.')
        else
            errorbar(J-1/4,mds(j),2*sds(j),'.b','markersize',28,'HandleVisibility','off')
        end
        hold on

        if J == 1
            errorbar(J+1/4,mdr(j),2*sdr(j),'.r','markersize',28,'DisplayName','R breakpoint difference mean \pm 2x s.e.')
        else
            errorbar(J+1/4,mdr(j),2*sdr(j),'.r','markersize',28,'HandleVisibility','off')
        end
        hold on
    end
    
    
    plot([1/2 length(Drugs)+1/2],[0 0],'-k','HandleVisibility','off');
    axis tight
    legend('location','southeast');
    set(gca,'Xtick',1:length(Drugs))
    set(gca,'Xticklabels',Drugs(srtmdsj))
    text(1.5,-4,'+CLSI bias')
    text(1.5,4,'+EUCAST bias')
    ylim([-5 5])
    set(gca,'Ytick',-5:5)
    set(gca,'Yticklabels',[5 4 3 2 1 0 1 2 3 4 5])
    ylabel('abs MIC difference (log2 \mug/mL)')
end
