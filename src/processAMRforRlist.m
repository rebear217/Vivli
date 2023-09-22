function processAMRforRlist(AMRforRMatches)
    close all

    P = globalParameterValues();
    lDosages = P.lDosages;

    ESRstr = AMRforRMatches.EUCASTbreakpointsSR;
    CSRstr = AMRforRMatches.CLSIbreakpointsSR;

    eucastS = zeros(1,length(ESRstr)-1);
    eucastR = zeros(1,length(ESRstr)-1);
    clsiS = zeros(1,length(ESRstr)-1);
    clsiR = zeros(1,length(ESRstr)-1);
    genus = zeros(1,length(ESRstr)-1);
    genusCount = 1;
    oldGenus = 'Acinetobacter';
    genusList = {};
    for J = 2:length(ESRstr)
        j = J-1;
        thisGenus = split(char(AMRforRMatches.Org(J)),' ');
        thisGenus = thisGenus{1};
        if strcmp(thisGenus(end),',')
            thisGenus = thisGenus(1:end-1);
        end

        if strcmp(oldGenus,thisGenus)
            genus(j) = genusCount;
        else
            genusCount = genusCount + 1;
            genus(j) = genusCount;
            genusList = [genusList(:)',oldGenus];
        end
        oldGenus = thisGenus;

        ESRsplit = split(ESRstr{J},"-");
        eucastS(j) = log2(str2num(ESRsplit{1}));
        eucastR(j) = log2(str2num(ESRsplit{2}));
        
        CSRsplit = split(CSRstr{J},"-");
        clsiS(j) = log2(str2num(CSRsplit{1}));
        clsiR(j) = log2(str2num(CSRsplit{2}));
    end
    genusList = [genusList(:)',thisGenus];

    %for j = 1:length(genusList)
    %    thisGenus = genusList{j};
    %    if strcmp(thisGenus(end),',')
    %        thisGenus = thisGenus(1:end-1);
    %    end
    %    genusList{j} = thisGenus;
    %end

    figure(1)
    shortLegend = 0;
    basicPlot(clsiS,eucastS,'S','.b');
    p = basicPlot(clsiR,eucastR,'R','.r');
    delete(p);

    figure(2)
    basicPlot([clsiS clsiR],[eucastS eucastR],'S and R');

    figure(3)
    set(3,'pos',[74         408        2487         929])

    ug = unique(genus);
    shortLegend = 1;
    for k = 1:length(ug)
        Fk = find(genus == k);
        cS = clsiS(Fk);
        cR = clsiR(Fk);
        eS = eucastS(Fk);
        eR = eucastR(Fk);
        subplot(3,9,k);

        try
            p = basicPlot(cS,eS,'S','.b');
            delete(p);
        end
        try
            basicPlot(cR,eR,'R','.r');
        catch
            p=plot([-10 10],[-10 10],'-k','DisplayName','equality');
        end

        title(genusList{k})
    end

    function p=basicPlot(X,Y,SR,marker)
        if nargin < 4
            marker = '.b';
        end
        X = X(~isnan(X));
        Y = Y(~isnan(X));
        X = X(~isnan(Y));
        Y = Y(~isnan(Y));
    
        Y = Y(X>-9);
        X = X(X>-9);
        X = X(Y>-9);
        Y = Y(Y>-9);

        Y = Y(X<9);
        X = X(X<9);
        X = X(Y<9);
        Y = Y(Y<9);

        rho = corr(X(:),Y(:));

        if shortLegend
            plot(X,Y,marker,'markersize',16,'DisplayName',[SR,', \rho \approx ',num2str(rho,3)]);
            xlabel('CLSI breakpoint')
            ylabel('EUCAST breakpoint')
            hold on
            p=plot([-10 10],[-10 10],'-k','DisplayName','equality');
            scale = 1.5;
            set(gca,'Xtick',-9:9)
            set(gca,'XtickLabel',-9:9)
            set(gca,'Ytick',-9:9)
            set(gca,'YtickLabel',-9:9)
        else
            plot(X,Y,marker,'markersize',16,'DisplayName',[SR,' breakpoints Pearson correlation \rho \approx ',num2str(rho,3)]);
            xlabel('CLSI breakpoint (log2 \mug/mL)')
            ylabel('EUCAST breakpoint  (log2 \mug/mL)')
            hold on
            p=plot([-10 10],[-10 10],'-k','DisplayName','EU-US equality');
            scale = 0.5;
            set(gca,'Xtick',lDosages)
            set(gca,'XtickLabel',round(lDosages,2))
            set(gca,'Ytick',lDosages)
            set(gca,'YtickLabel',round(lDosages,2))
        end
        xlim([-10 10])
        ylim([-10 10])
        legend('location','southwest');
        sizePlot(X,Y,scale,marker);
        grid on
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