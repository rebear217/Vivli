function [synthdata,bestfit,R2] = ECOFFexamples()

P = globalParameterValues;
close all
ms = 34;

for j = 1:2
    switch j
        case 1
            [synthdata,bins] = syntheticData([200 10 30],[-6 1 4],[1 1 1]);
        case 2
            [synthdata,bins] = syntheticData([60 1 230],[-6 1 4],[1 1 1]);
    end

    [bestfit,bestbin,R2,~,ECOFF,~,allFits] = findECOFF(synthdata,bins,P.ECOFFpValue);    

    figure(j)
    set(j,'pos',[39   224   813   473])

    plot(bins,synthdata/sum(synthdata),'.--k','markersize',ms,'linewidth',1,'DisplayName','MIC histogram')
    hold on
    p=plot(bins,bestfit.feval(bins),'-k','linewidth',4,'DisplayName',['best datasize is ',num2str(bestbin)]);
    %title('best R^2 fit')
    ylabel('frequency')
    xlabel('MIC (log2 \mug/mL)')
    plot(ECOFF,0,'or','markersize',30,'DisplayName','ECOFF')
    axis tight

    for k = 2:length(allFits)
        s = (k-1) / (length(allFits)-1);
        colr = [1-s 0 s];
        muguess = allFits{k}.Coefficients.Estimate(1);
        sigmaguess = allFits{k}.Coefficients.Estimate(2);    
        plot(bins,gaussian(bins,muguess,sigmaguess),'-',...
            'color',colr,'LineWidth',2,'DisplayName',[num2str(length(allFits{k}.Fitted)),' datapoints'])
        %subplot(1,2,2)
        %plot(length(allFits{k}.Fitted),allFits{k}.Rsquared.Ordinary,'o',...
        %    'color',colr,'markersize',20)
        %subplot(1,2,1)
    end
    legend('location','northeast')
    
    %{
    if j == 1
        ax = axes('Position',[0.45 0.7 0.25 0.25]);
    else
        ax = axes('Position',[0.2 0.6 0.35 0.35]);
    end
    %subplot(1,2,2)
    plot(R2,'.-b','markersize',ms,'linewidth',1,'DisplayName','EUCAST datafit criterion (1 is best)');
    [mj,m] = max(R2);
    hold on
    plot(m,mj,'or','markersize',ms,'DisplayName','optimal datasize')
    ylabel('R^2')
    xlabel('datasize attempted')
    axis tight
    yL = ylim;
    ylim([yL(1) yL(2)*1.1]);
    %title('maximise goodness of normal fit')
    legend()   
    %}
    
end

