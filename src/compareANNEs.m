function compareANNEs(CnewHistosANNdecisionData,EnewHistosANNdecisionData)

    bins = -9:0.125:9;

    figure(1)
    X = CnewHistosANNdecisionData.table.ANNbpS;
    Y = EnewHistosANNdecisionData.table.ANNbpS;
    s = skewness(X-Y);
    sp = skekurtest(X-Y,1);

    plot(X,Y,'.b','DisplayName',['S breakpoint predictions (\rho \approx ',num2str(corr(X,Y),3),')'],'markersize',14)
    
    figure(2)
    histogram(X-Y,'BinEdges',bins,'Normalization','probability','DisplayName',['S (skew ',num2str(s,2),', p \approx ',num2str(sp,2),')'])

    figure(1)
    hold on
    X = CnewHistosANNdecisionData.table.ANNbpR;
    Y = EnewHistosANNdecisionData.table.ANNbpR;
    plot(X,Y,'.r','DisplayName',['R breakpoint predictions (\rho \approx ',num2str(corr(X,Y),3),')'],'markersize',14)

    figure(2)
    hold on
    s = skewness(X-Y);
    sp = skekurtest(X-Y,1);
    histogram(X-Y,'BinEdges',bins,'Normalization','probability','DisplayName',['R (skew ',num2str(s,2),', p \approx ',num2str(sp,2),')'])
    xlim([-3 3])
    ylabel('frequency')
    xlabel('ANNE breakpoint differences (CLSI-EUCAST;log2 \mug/mL)')
    legend();

    figure(1)
    x = [-9 9];
    plot(x,x,'-k','DisplayName','CLSI-EUCAST equality')
    xlabel('CLSI ANNE decisions (log2 \mug/mL)')
    ylabel('EUCAST ANNE decisions (log2 \mug/mL)')
    xlim([-4 9])
    ylim([-4 9])
    
    legend('location','northwest')
    grid on

    figure(3)
    X = CnewHistosANNdecisionData.table.ANNbpSIQRhigh - CnewHistosANNdecisionData.table.ANNbpSIQRlow;
    Y = EnewHistosANNdecisionData.table.ANNbpSIQRhigh - EnewHistosANNdecisionData.table.ANNbpSIQRlow;
    s = skewness(Y-X);
    plot(X,Y,'.b','DisplayName',['S IQR'],'markersize',14)

    hold on

    X = CnewHistosANNdecisionData.table.ANNbpRIQRhigh - CnewHistosANNdecisionData.table.ANNbpRIQRlow;
    Y = EnewHistosANNdecisionData.table.ANNbpRIQRhigh - EnewHistosANNdecisionData.table.ANNbpRIQRlow;
    s = skewness(Y-X);
    plot(X,Y,'.r','DisplayName',['R IQR'],'markersize',14)    

    plot(x,x,'-k','DisplayName','CLSI-EUCAST equality')
    xlabel('CLSI ANNE decision IQRs (log2 \mug/mL)')
    ylabel('EUCAST ANNE decision IQRs (log2 \mug/mL)')
    xlim([0 9])
    ylim([0 9])
    legend('location','northwest')

end