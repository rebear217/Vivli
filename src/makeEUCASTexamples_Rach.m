function [synthdata,bestfit,R2] = makeEUCASTexamples_Rach()

P = globalParameterValues;
close all

for j = 1:2
    switch j
        case 1
            [synthdata,bins] = syntheticData([200 10 30],[-6 1 4],[1 1 1]);
        case 2
            [synthdata,bins] = syntheticData([60 1 230],[-6 1 4],[1 1 1]);
    end

    [bestfit,bestbin,R2,wholeRMS,ECOFF] = findECOFF(synthdata,bins,P.ECOFFpValue);
    

    figure(j)
    set(j,'pos',[62   674   896   663])

    subplot(2,2,1)
    plot(bins,synthdata,'.k','markersize',20,'DisplayName','synth data');
    ylabel('histogram')
    xlabel('log2 MIC dose')
    axis tight
    legend()
    title('a synthetic histogram')

    subplot(2,2,2)
    plot(R2,'.-b','markersize',20,'linewidth',1,'DisplayName','EUCAST datafit criterion (1 is best)');
    [mj,m] = max(R2);
    hold on
    plot(m,mj,'or','markersize',20,'DisplayName','optimum')
    ylabel('R^2')
    xlabel('datasize attempted')
    axis tight
    yL = ylim;
    ylim([yL(1) 1]);
    title('maximise goodness of normal fit')
    legend()
    
    subplot(2,2,3)
    plot(wholeRMS,'.r','markersize',20);
    ylabel('entire RMS error')
    xlabel('datasize attempted')
    axis tight
    title('error from entire histogram')

    subplot(2,2,4)
    plot(bins,synthdata/sum(synthdata),'.k','markersize',20)
    hold on
    p=plot(bins,bestfit.feval(bins),'-k','linewidth',1);
    title('best R^2 fit')
    ylabel('frequency')
    xlabel('log2 MIC dose')
    legend(p,['best datasize is ',num2str(bestbin)])
    plot(ECOFF,0,'or','markersize',20,'DisplayName','ECOFF')
    axis tight

    export_fig(['./figures/EUCASTtest',num2str(j),'.pdf'])

end

