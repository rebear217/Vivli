function ANNanalyseOneHistogram(drug,bug,ANNdecisionData,CEflag,figNum)
    
    if nargin < 5
        figNum = 1;
    end

    P = globalParameterValues();
    bins = min(P.lDosages):0.33:max(P.lDosages);

    %analyseOneHistogram(allEUCASTHistograms,drug,bug);

    Fd = find(strcmpi(ANNdecisionData.table.drugs,drug));
    Fb = find(strcmpi(ANNdecisionData.table.bugs,bug));
    F = intersect(Fd,Fb);
    disp(['Found at position ',num2str(F)]);
    
    if ~isempty(F)
        
        if CEflag
            bpS = ANNdecisionData.table.CLSIbpS(F);
            bpR = ANNdecisionData.table.CLSIbpR(F);
        else
            bpS = ANNdecisionData.table.EUCASTbpS(F);
            bpR = ANNdecisionData.table.EUCASTbpR(F);
        end
    
        ANNbpS = ANNdecisionData.table.ANNbpS(F);
        ANNbpR = ANNdecisionData.table.ANNbpR(F);
        
        ANNRlowIQR = ANNdecisionData.table.ANNbpRIQRlow(F);
        ANNRhighIQR = ANNdecisionData.table.ANNbpRIQRhigh(F);
        
        ANNSlowIQR = ANNdecisionData.table.ANNbpSIQRlow(F);
        ANNShighIQR = ANNdecisionData.table.ANNbpSIQRhigh(F);
    
        ANNRdecisions = ANNdecisionData.ensembleRDecisions(F,:);
        ANNSdecisions = ANNdecisionData.ensembleSDecisions(F,:);
    
        ANNhistogram = ANNdecisionData.histograms(F,:);
    
        figure(figNum)
        set(figNum,'pos',[463   770   838   458])
    
        plot(P.lDosages,ANNhistogram,'.-k','MarkerSize',34,'DisplayName',['MIC frequency data (',bug,' & ',drug,')'])
        hold on
        histogram(ANNRdecisions,bins,'Normalization','probability','FaceColor','r','DisplayName','ANN ensemble Rs')
        histogram(ANNSdecisions,bins,'Normalization','probability','FaceColor','b','DisplayName','ANN ensemble Ss')

        axis tight;
        yL = ylim;
        M = 1.3*yL(2);
        ylim([0 1.4*yL(2)])
    
        plot([ANNbpS ANNbpS],[0 0.95*M],'-b','linewidth',1,'DisplayName','ANN mean S')
        plot([ANNbpR ANNbpR],[0 0.95*M],'-r','linewidth',1,'DisplayName','ANN mean R')
        plot([ANNRlowIQR ANNRhighIQR],[1 1]*0.9*M,'-r','linewidth',4,'DisplayName',['ANN ensemble R IQR (',num2str(round(2^ANNRlowIQR,1)),',',num2str(round(2^ANNRhighIQR,1)),')\mug/mL']);
        plot([ANNSlowIQR ANNShighIQR],[1 1]*0.9*M,'-b','linewidth',4,'DisplayName',['ANN ensemble S IQR (',num2str(round(2^ANNSlowIQR,1)),',',num2str(round(2^ANNShighIQR,1)),')\mug/mL']);

        if CEflag
            plot([bpS bpS],[0 0.95*M],'--b','linewidth',3,'DisplayName','CLSI breakpoint S')
            plot([bpR bpR],[0 0.95*M],'--r','linewidth',3,'DisplayName','CLSI breakpoint R')
        else
            plot([bpS bpS],[0 0.95*M],'--b','linewidth',3,'DisplayName','EUCAST breakpoint S')
            plot([bpR bpR],[0 0.95*M],'--r','linewidth',3,'DisplayName','EUCAST breakpoint R')
        end
    
        xlabel('MIC (log2 \mug/mL)')
        ylabel('frequency')
        legend('location','northwest')
        legend('boxoff')
    else
        M = length(ANNdecisionData.table.drugs);
        disp(['No matches found for that drug and bug where the ANN decisions dataset has ',num2str(M),' many elements.'])
    end
end

