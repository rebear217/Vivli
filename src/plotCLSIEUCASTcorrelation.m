function plotCLSIEUCASTcorrelation(CLSIEUCASTMatches)
    P = globalParameterValues();
    lDosages = P.lDosages;

    x = CLSIEUCASTMatches(:,6);
    y = CLSIEUCASTMatches(:,8);
    X = zeros(1,length(x));
    Y = zeros(1,length(y));
    
    for j = 1:length(x)
        try
            X(j) = x{j};
            Y(j) = y{j};        
        catch
            X(j) = NaN;
            Y(j) = NaN;        
        end
    end
    
    %%
    
    XX = X(~isnan(X));
    YY = Y(~isnan(X));
    X = XX(~isnan(YY));
    Y = YY(~isnan(YY));
    
    figure(1)
    rho = corr(X(:),Y(:));
    plot(X,Y,'.b','markersize',8,'DisplayName',['Pearson correlation \rho \approx ',num2str(rho)]);
    hold on
    plot([-10 10],[-10 10],'-k','DisplayName','equality')
    xlim([-10 10])
    ylim([-10 10])
    set(gca,'Xtick',lDosages)
    set(gca,'XtickLabel',round(lDosages,2))
    set(gca,'Ytick',lDosages)
    set(gca,'YtickLabel',round(lDosages,2))
    xlabel('CLSI breakpoint S')
    ylabel('EUCAST breakpoint S')
    legend('location','northwest');
    sizePlot(X,Y);

    %%
    
    x = CLSIEUCASTMatches(:,7);
    y = CLSIEUCASTMatches(:,9);
    X = zeros(1,length(x));
    Y = zeros(1,length(y));
    
    for j = 1:length(x)
        try
            X(j) = x{j};
            Y(j) = y{j};        
        catch
            X(j) = NaN;
            Y(j) = NaN;        
        end
    end
    
    %%
    
    XX = X(~isnan(X));
    YY = Y(~isnan(X));
    X = XX(~isnan(YY));
    Y = YY(~isnan(YY));
    
    figure(2)
    rho = corr(X(:),Y(:));
    plot(X,Y,'.b','markersize',8,'DisplayName',['Pearson correlation \rho \approx ',num2str(rho)]);
    hold on
    plot([-10 10],[-10 10],'-k','DisplayName','equality')
    xlim([-10 10])
    ylim([-10 10])
    set(gca,'Xtick',lDosages)
    set(gca,'XtickLabel',round(lDosages,2))
    set(gca,'Ytick',lDosages)
    set(gca,'YtickLabel',round(lDosages,2))
    xlabel('CLSI breakpoint R')
    ylabel('EUCAST breakpoint R')
    legend('location','northwest');
    sizePlot(X,Y);

    function sizePlot(X,Y,scale)
        if nargin < 3
            scale = 3;
        end
        UX = unique(X);
        UY = unique(Y);
        for j = 1:length(UX)
            for i = 1:length(UY)
                S = sum( (X == UX(j)) & (Y == UY(i)) );
                if S > 0
                    plot(UX(j),UY(i),'.b','markersize',scale*S,'HandleVisibility','off');
                end
            end
        end
    end

end
