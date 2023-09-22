function illustrateNewHistoDecisions(newHistosANNdecisionData,CEflag)

    P = globalParameterValues();
    lDosages = P.lDosages;
    Nclusters = 25;

    figure(1)
    set(1,'pos',[138         268        1241         570])
    idx = weightOrdering(newHistosANNdecisionData.histograms);
    imagesc(flipud(newHistosANNdecisionData.histograms(idx,:)'));
    hold on

    title('Heatmap of observed frequencies')
    colormap(1-bone)
    ylabel('MIC dose (log2 \mug/mL)')
    set(gca,'Ytick',1:19)
    set(gca,'Yticklabel',fliplr(round(lDosages,2)))
    xlabel('spectrally ordered EUCAST histograms with ANN breakpoints')
    colorbar
    drawnow

    [N,~] = size(newHistosANNdecisionData.histograms);
    colr = [0.5 0.5 1];
    for i = 1:N
        if ~CEflag
            check = newHistosANNdecisionData.table.EUCASTbpS(idx(i));
        else
            check = 0;
        end

        if isnan(check) || (check > -9)
            L = 10-newHistosANNdecisionData.table.ANNbpS(idx(i));
            R = 10-newHistosANNdecisionData.table.ANNbpR(idx(i));
            p1 = plot([i i],[L R],'-','linewidth',1,'color',colr);
    
            if CEflag
                L = 10-newHistosANNdecisionData.table.CLSIbpS(idx(i));
                R = 10-newHistosANNdecisionData.table.CLSIbpR(idx(i));
            else
                L = 10-newHistosANNdecisionData.table.EUCASTbpS(idx(i));
                R = 10-newHistosANNdecisionData.table.EUCASTbpR(idx(i));
            end
            p2 = plot([i i],[L R],'-r','linewidth',1);
        end
    end

    p3 = plot(-20,20,'-k');
    if CEflag
        legend([p1,p2,p3],{'ANN breakpoints','CLSI breakpoints','histograms'},'location','southeast');
    else
        legend([p1,p2,p3],{'ANN breakpoints','EUCAST breakpoints','histograms'},'location','southeast');
    end

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
    
    Xrange = [-8 10];
    for j = 1:2
        figure(j+1)
        if j == 1
            if CEflag
                X = newHistosANNdecisionData.table.CLSIbpS;
                xStr = 'CLSI S breakpoints (log2 \mug/mL)';
            else
                X = newHistosANNdecisionData.table.EUCASTbpS;
                xStr = 'EUCAST S breakpoints (log2 \mug/mL)';
            end
            Y = newHistosANNdecisionData.table.ANNbpS;
            yStr = 'ANN S breakpoints (log2 \mug/mL)';
        else
            if CEflag
                X = newHistosANNdecisionData.table.CLSIbpR;
                xStr = 'CLSI R breakpoints (log2 \mug/mL)';
            else
                X = newHistosANNdecisionData.table.EUCASTbpR;
                xStr = 'EUCAST R breakpoints (log2 \mug/mL)';
            end
            Y = newHistosANNdecisionData.table.ANNbpR;
            yStr = 'ANN R breakpoints (log2 \mug/mL)';
        end

        Y = Y(~isnan(X));
        X = X(~isnan(X));

        X = X(~isnan(Y));
        Y = Y(~isnan(Y));

        if ~CEflag
            %eucase has some weird -ve data in there, so don't include it:
            Y = Y(X > -9);
            X = X(X > -9);

            X = X(Y > -9);
            Y = Y(Y > -9);            
        end

        lf = fitlm(X,Y);
        plot(X,Y,'.k','markersize',20);
        hold on
        p2 = plot(Xrange,Xrange,'-k');
        p3 = plot(Xrange,lf.feval(Xrange),'--k','linewidth',1);

        xlabel(xStr)
        ylabel(yStr)
        set(gca,'Xticklabel',round(lDosages,2))
        set(gca,'Yticklabel',round(lDosages,2))
        xlim(Xrange)
        ylim(Xrange)
        if CEflag
            legend([p2,p3],{'ANN-CLSI equality',['linear regression (\rho\approx ',num2str(corr(X,Y),3),')']})
        else
            legend([p2,p3],{'ANN-EUCAST equality',['linear regression (\rho\approx ',num2str(corr(X,Y),3),')']})
        end
        axis tight
        grid on
        
    end

end