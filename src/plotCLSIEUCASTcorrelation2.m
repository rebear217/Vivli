function plotCLSIEUCASTcorrelation2(bothBreakpointSets)

    close all
    drawArrow = @(x,y,varargin) quiver( x(1),y(1),x(2)-x(1),y(2)-y(1),0,varargin{:});

    P = globalParameterValues();
    lDosages = P.lDosages;

    figure(1)
    XS = bothBreakpointSets.S_clsi;
    YS = bothBreakpointSets.S_eucast;
    XS = XS(YS > -9);
    YS = YS(YS > -9);
    
    basicPlot(XS,YS,'S','.b');

    XR = bothBreakpointSets.R_clsi;
    YR = bothBreakpointSets.R_eucast;
    XR = XR(YR > -9);
    YR = YR(YR > -9);

    p=basicPlot(XR,YR,'R','.r');
    delete(p);
    grid on

    figure(3)
    XB = [XS ; XR];
    YB = [YS ; YR];
    basicPlot(XB,YB,'S and R');
    grid on

    figure(4);
    skS = skewness(XS-YS);
    skR = skewness(XR-YR);
    skekurtestSp = skekurtest(XS-YS,1);
    skekurtestRp = skekurtest(XR-YR,1);

    strS = [num2str(mean(XS-YS),3),'\pm',num2str(ste(XS-YS),3)];
    strR = [num2str(mean(XR-YR),3),'\pm',num2str(ste(XR-YR),3)];

    histogram(XS-YS,'DisplayName',['S (mean ',strS,'s.e., skew ',num2str(skS,2),', p\approx',num2str(skekurtestSp,2),')']);
    hold on
    histogram(XR-YR,'DisplayName',['R (mean ',strR,'s.e., skew ',num2str(skR,2),', p\approx',num2str(skekurtestRp,2),')']);
    axis tight
    xlabel('breakpoint differences (log2 \mug/mL)')
    ylabel('#observations')
    legend('location','northwest')
    %drawArrow([0,8],[80,80],'-k','linewidth',1,'DisplayName','US > EU breakpoints')

    function p=basicPlot(X,Y,SR,marker)
        if nargin < 4
            marker = '.b';
        end
        X = X(~isnan(X));
        Y = Y(~isnan(X));
        X = X(~isnan(Y));
        Y = Y(~isnan(Y));

        %Y = Y(X>-9);
        %X = X(X>-9);
        %X = X(Y>-9);
        %Y = Y(Y>-9);

        rho = corr(X,Y);
        plot(X,Y,marker,'markersize',16,'DisplayName',[SR,' breakpoints Pearson correlation \rho \approx ',num2str(rho)]);
        hold on
        p=plot([-10 10],[-10 10],'-k','DisplayName','EU-US equality');
        xlim([-10 10])
        ylim([-10 10])
        set(gca,'Xtick',lDosages)
        set(gca,'XtickLabel',round(lDosages,2))
        set(gca,'Ytick',lDosages)
        set(gca,'YtickLabel',round(lDosages,2))
        xlabel('CLSI breakpoint (log2 \mug/mL)')
        ylabel('EUCAST breakpoint  (log2 \mug/mL)')
        legend('location','southwest');
        sizePlot(X,Y,3,marker);
    end

    function sizePlot(X,Y,scale,marker)
        UX = unique(X);
        UY = unique(Y);
        for j = 1:length(UX)
            for i = 1:length(UY)
                S = sum( (X == UX(j)) & (Y == UY(i)) );
                if S > 0
                    plot(UX(j),UY(i),marker,'markersize',scale*S,'HandleVisibility','off');
                end
            end
        end
    end
end