function analyseOneHistogram(allEUCASTHistograms,drug,bug,plotOFF)
    
    close all

    Params = globalParameterValues();
    foptions = optimoptions('fsolve','Display','none');
    
    if nargin < 4
        plotOFF = 0;
    end

    [N,M] = size(allEUCASTHistograms);
    found = zeros(N,1);
    for j = 1:N
        thisDrug = allEUCASTHistograms{j,1};
        if length(thisDrug) > length(drug)
            thisDrug = thisDrug(1:length(drug));
        end
        thisBug = allEUCASTHistograms{j,2};
        if length(thisBug) > length(bug)
            thisBug = thisBug(1:length(bug));
        end
        if strcmpi(thisBug,bug) && strcmpi(thisDrug,drug)
            found(j) = 1;
        end
    end

    F = find(found);
    if isempty(F)
        disp('No drug-bug pairs match your text')
    else
        if length(F) > 1
            disp('------------------------------------------------------------')
            disp('Found more than one matching drug-bug pair, using the first:')
            disp('------------------------------------------------------------')
            for j = 1:length(F)
                disp([allEUCASTHistograms{F(j),1},' & ',allEUCASTHistograms{F(j),2}])
            end
            F = F(1);
        end

        MICs = cell2mat(allEUCASTHistograms(F,3:M));

        ECOFFpValue = Params.ECOFFpValue;
        noiseLevel = Params.noiseLevel;        
        Ngmmodels = Params.Ngmmodels;
        options = statset('MaxIter',2000);
        GMModel = cell(Ngmmodels,1);
        lDosages = Params.lDosages;
        lDosagesString = Params.lDosagesString;
        finelDosages = min(lDosages):0.01:max(lDosages);
        fineDosages = 2.^finelDosages;

        disp(' ')
        disp('-- ECOFF Details ----------------------------------------------------------')
        disp(' ')

        try
            [Xbestfit,Xbestbin,XR2,XwholeRMS,XECOFF,XECOFFCorr] = findECOFF(MICs,lDosages,ECOFFpValue);
            disp(['Found (T)ECOFF at ',num2str(XECOFF),' log2 ug/mL']);
            disp(['Histogram correlation below ECOFF is ',num2str(XECOFFCorr),' where 1 is best, -1 is worst']);
            disp(['Best R^2 is ',num2str(XR2(Xbestbin)),' where 1 is best, 0 is worst']);
            disp(' ')
            disp('Optimal ECOFF Gaussian fit details are ...')
            disp(Xbestfit)
        catch
            disp('Algorithm cannot determine ECOFF')
            XECOFF = NaN;
            XECOFFCorr = NaN;
        end

        disp(' ')
        disp('-- GMM Details ------------------------------------------------------------')
        disp(' ')

        [newdist,allnewPatientMICs] = generatePatientsMICsNOISE(MICs,lDosages,noiseLevel);
        snd = sum(newdist);

        AIC = NaN(1,Ngmmodels);
        for k = 1:Ngmmodels
            try
                %GMModel{k} = fitgmdist(allnewPatientMICs',k,'Options',options,'RegularizationValue',0.1);
                GMModel{k} = fitgmdist(allnewPatientMICs',k,'Options',options);
                AIC(k) = GMModel{k}.AIC;
            catch
                AIC(k) = NaN;
            end
        end
        [m,I] = min(AIC);
    
        gmm = GMModel{I};
        f = @(x)gmm.pdf(x);
        Post = @(x)posterior(gmm,x);
        
        Rsquared = 1 - sum((snd*f(lDosages') - MICs').^2) / sum((MICs - mean(MICs)).^2);

        figure(1)
        set(1,'pos',[68   375   598   322])
        plot(lDosages,MICs,'.k','markersize',20,'DisplayName',...
            [bug,' & ',drug,' MIC histogram data']);
        hold on
        plot(finelDosages',snd*f(finelDosages'),'-b','linewidth',2,...
            'DisplayName',['GMM (',num2str(I),' clusters)']);
        %title([allEUCASTHistograms{F,1},' & ',allEUCASTHistograms{F,2}])
        legend boxoff
        legend('location','northwest')
        xlabel('MIC (log2 \mug/mL)')
        ylabel('no. of cases')
        set(gca,'Xtick',round(lDosages,2))
        set(gca,'Xticklabels',round(lDosages,2))
        axis tight
        yL = ylim;
        ylim([yL(1) 1.2*yL(2)])

        [mu1,muI] = min(gmm.mu);
        sigma1 = gmm.Sigma(muI);
        gmmCP = gmm.ComponentProportion(muI);
        g = gmdistribution(mu1,sigma1);
        G = @(x)gmmCP*(g.pdf(x));
    
        [mu2,muI2] = max(gmm.mu);
        sigma2 = gmm.Sigma(muI2);
        gmmCP2 = gmm.ComponentProportion(muI2);
        g2 = gmdistribution(mu2,sigma2);
        G2 = @(x)gmmCP2*(g2.pdf(x));

        if I > 2
            [Muu,MI] = mink(gmm.mu,3);
            sigmaInt = gmm.Sigma(MI(2));
            muInt = gmm.mu(MI(2));
            gmmCPInt = gmm.ComponentProportion(MI(2));
            gInt = gmdistribution(muInt,sigmaInt);

            GInt = @(x)gmmCPInt*(gInt.pdf(x));        
            SIboundary = fsolve(@(x)(log10(G(x))-log10(GInt(x))),muInt,foptions);
            IRboundary = fsolve(@(x)(log10(G2(x))-log10(GInt(x))),(muInt+mu2)/2,foptions);
            
            text(Muu(1),0.1*yL(2),'S','fontsize',22,'FontName','times')
            text(Muu(2),0.1*yL(2),'I','fontsize',22,'FontName','times')
            text(Muu(3),0.1*yL(2),'R','fontsize',22,'FontName','times')
        end

        figure(2)
        set(2,'pos',[68   375   598   322])

        if gmm.Converged
            bLabel = [bug,' & ',drug,' GMM (',num2str(I),' clusters)'];
        else
            bLabel = [bug,' & ',drug,' (',num2str(I),' clusters w/o converging)'];
        end
        plot(finelDosages',f(finelDosages'),'-b','linewidth',2,'DisplayName',bLabel);
        hold on
        plot(finelDosages',G(finelDosages'),'-k','linewidth',2,'DisplayName','S');
        axis tight
        yL = ylim;

        if I > 2
            plot(finelDosages',GInt(finelDosages'),'-r','linewidth',1,'DisplayName','I');
        end
        if I > 1
            plot(finelDosages',G2(finelDosages'),'-r','linewidth',2,'DisplayName','R');
        end
        if I > 2
            plot([SIboundary SIboundary],[yL(1) yL(2)],'--b','DisplayName','SI boundary','linewidth',2)    
            plot([IRboundary IRboundary],[yL(1) yL(2)],'--k','DisplayName','IR boundary','linewidth',2)
        end        
        if I > 1
            SRboundary = fsolve(@(z)(log2(G(z))-log2(G2(z))),(mu1+mu2)/2,foptions);        
            plot([SRboundary SRboundary],[yL(1) yL(2)],'--r','DisplayName','SR boundary','linewidth',2)
        end
        if I == 2
            checkSRb = solveForWTBoundary(Post,1,2,SRboundary);
            disp('--------------------------------------------')
            disp(['SR boundary (equal distributions) : ',num2str(SRboundary)])
            disp(['Point with equal posteriors : ',num2str(checkSRb)])
            disp('--------------------------------------------')
        end
        if I == 3
            checkSRb = [0,0,0];
            checkSRb(1) = solveForWTBoundary(Post,MI(1),MI(2),SIboundary);
            checkSRb(2) = solveForWTBoundary(Post,MI(2),MI(3),IRboundary);
            checkSRb(3) = solveForWTBoundary(Post,MI(1),MI(3),SRboundary);
            disp('--------------------------------------------')
            disp(['Points with equal distributions : ',num2str([SIboundary IRboundary SRboundary])])
            disp(['Points with equal posteriors : ',num2str(checkSRb)])
            disp('--------------------------------------------')
        end

        if ~isnan(XECOFF)
            plot(XECOFF,0,'.k','linewidth',1,'markersize',34,'DisplayName','est. ECOFF');
        end
        %title([allEUCASTHistograms{F,1},' & ',allEUCASTHistograms{F,2}])

        %plot(lDosages,MICs/sum(MICs),'.k','markersize',32,'DisplayName','MIC histogram data','color',[1 1 1]*0.6);

        legend boxoff
        legend('location','northwest')
        axis tight
        yL = ylim;
        ylim([yL(1) 1.15*yL(2)])
        xlabel('MIC (log2 \mug/mL)')
        ylabel('likelihood')
        set(gca,'Xtick',round(lDosages,2))
        set(gca,'Xticklabels',round(lDosages,2))

        RLmatrix = exp(-abs((AIC - AIC')));
        Lmatrix = log10(RLmatrix - eye(Ngmmodels));
        matrixMax = max(max(Lmatrix));

        disp(['Did GMM converge?: ',num2str(gmm.Converged)])
        disp(['Optimal Gauss mixture number: ',num2str(I),' bumps'])
        disp(['Gauss mixture R^2: ',num2str(Rsquared),' where 1 is best, 0 is worst'])
        disp(['Next best GMM log10 likelihood: ',num2str(matrixMax)])
        disp(' ');
        disp(gmm);
        

    end

end

