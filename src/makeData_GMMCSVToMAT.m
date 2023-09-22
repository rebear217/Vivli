clear all
close all
clc
warning('off')

%%

Params = globalParameterValues();

PDFsave = 0;
%ignoreList = {'Ceftazidime-clavulanate.csv'};
ignoreList = [];

breakpoints = importCLSIbreakpoint('./CLSI/breakpoints.csv');
ECOFFpValue = Params.ECOFFpValue;

disp(['Imported ',num2str(length(breakpoints.S)),' breakpoints'])

outputTitles = {'antibiotic','pathogen','optimal cluster number','GMM converged','R squared',...
    'S mean','S sigma','I mean','I sigma','R mean','R sigma',...
    'SI boundary','SR boundary','IR boundary','max log10 relative likelihood',...
    '(T)ECOFF','(T)ECOFF CI min','(T)ECOFF CI max','CLSI STIC S','CLSI STIC R',...
    '(X)ECOFF','(X)ECOFF correlation','EUCAST S bp','EUCAST R bp','GMMweightS','GMMweightI','GMMweightR'};
outputUnits = {'string','string','integer','Boolean','real',...
    'log2 ug/mL','(log2 ug/mL)^2','log2 ug/mL','(log2 ug/mL)^2','log2 ug/mL','(log2 ug/mL)^2',...
    'log2 ug/mL','log2 ug/mL','log2 ug/mL','real',...
    'log2 ug/mL','log2 ug/mL','log2 ug/mL','log2 ug/mL','log2 ug/mL',...
    'log2 ug/mL','real','log2 ug/mL','log2 ug/mL','real','real','real'};

%%

rootDir = './EUCASThistograms/';
d = dir([rootDir,'*.csv']);

analysisFile = d(1).name;

Dosages = importDosages([rootDir,analysisFile]);
lDosages = log2(Dosages);
finelDosages = min(lDosages):0.01:max(lDosages);
fineDosages = 2.^finelDosages;

totalFileNo = length(d);
entireSIRData = cell(3,length(outputTitles));
entireCount = 3;
consolidatedSIRData = cell(3,length(outputTitles));
consolidatedCount = 3;

foptions = optimoptions('fsolve','Display','none');

entirePDFs = {};

for fileJ = 1:totalFileNo
    
    close all

    analysisFile = d(fileJ).name;
    disp(analysisFile)

    loadError = 0;
    if ismember(analysisFile,ignoreList)
        loadError = 1;
    end
    try
        MICdata = importMICs([rootDir,analysisFile]);
        bugs = importBugNames([rootDir,analysisFile]);
        [TECOFF, ConfidenceInterval] = importECOFFs([rootDir,analysisFile]);
    catch
        disp('This has a non-standard CSV structure.')
        loadError = 1;
    end
    if isempty(bugs{1})
        loadError = 1;
        disp('This has a non-standard CSV structure (it loaded but with junk).')
    end

    if ~loadError
        TCnumbers = log2(TECOFFStoNumbers(TECOFF));
        [CIL,CIR] = CItoNumbers(ConfidenceInterval);
        CIL = log2(CIL);
        CIR = log2(CIR);
        
        %%

        thisDrugSp = getDrugFromFilename(analysisFile);
        disp(['    Check : ',analysisFile,' provides data for this drug ... ',thisDrugSp])

        %%
        
        [n,m] = size(MICdata);
        %P = round(sqrt(n));
        noiseLevel = 0;%no noise at all        
        Ngmmodels = Params.Ngmmodels;

        options = statset('MaxIter',2000);
        %replicate the GMM fits several times and pick the one with the best AIC:
        options.Replicates = Params.Replicates;
        GMModel = cell(Ngmmodels,1);
            
        outputCell = cell(n+2,length(outputTitles));
    
        [outputCell{1,:}] = deal(outputTitles{:});
        [outputCell{2,:}] = deal(outputUnits{:});
        
        [entireSIRData{1,:}] = deal(outputTitles{:});
        [entireSIRData{2,:}] = deal(outputUnits{:});
    
        [consolidatedSIRData{1,:}] = deal(outputTitles{:});
        [consolidatedSIRData{2,:}] = deal(outputUnits{:});
    
        %%
        
        for j = 1:n
            
            %this does not work: it creates far too many matches:
            %[CLSIS,CLISR] = getBreakpoint(bugs{j},thisDrugSp,breakpoints);
            [CLSIS,CLISR] = getBreakpointNEW(bugs{j},thisDrugSp,breakpoints);

            outputVector = {thisDrugSp,bugs{j},NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,CLSIS,CLISR,NaN,NaN,NaN,NaN,NaN,NaN,NaN};

            if PDFsave
                figure(j)
                set(j,'pos',[58         713        1403         384])
                subplot(1,3,1);
            end
    
            MICs = MICdata(j,:);
            try
                [Xbestfit,Xbestbin,XR2,XwholeRMS,XECOFF,XECOFFCorr] = findECOFF(MICs,lDosages,ECOFFpValue);
            catch
                XECOFF = NaN;
                XECOFFCorr = NaN;
            end
        
            outputVector{21} = XECOFF;
            outputVector{22} = XECOFFCorr;    
        
            [newdist,allnewPatientMICs] = generatePatientsMICsNOISE(MICdata(j,:),lDosages,noiseLevel);
            snd = sum(newdist);
        
            AIC = NaN(1,Ngmmodels);
            for k = 1:Ngmmodels
                try
                    GMModel{k} = fitgmdist(allnewPatientMICs',k,'Options',options,'RegularizationValue',0.1);
                    %GMModel{k} = fitgmdist(allnewPatientMICs',k,'Options',options);
                    AIC(k) = GMModel{k}.AIC;
                catch
                    AIC(k) = NaN;
                end
            end
        
            [m,I] = min(AIC);
            outputVector{3} = I;
        
            gmm = GMModel{I};
            f = @(x)gmm.pdf(x);
            
            if PDFsave
                plot(lDosages,MICs,'.k','markersize',20,'DisplayName',...
                    [bugs{j},' MICs']);
                hold on
                plot(finelDosages',snd*f(finelDosages'),'-b','linewidth',2,...
                    'DisplayName',['GMM (n=',num2str(I),')']);
            end
    
            Rsquared = 1 - sum((snd*f(lDosages') - MICs').^2) / sum((MICs - mean(MICs)).^2);
            outputVector{5} = Rsquared;
        
            if PDFsave && j == 1
                title(thisDrugSp)
            end
            if PDFsave
                legend();
                legend boxoff
                legend('location','northwest')
                xlabel('log2 MIC dose')
                ylabel('no. of cases')    
                axis tight
        
                subplot(1,3,2);
            end
        
            [muS,muIS] = min(gmm.mu);
            sigmaS = gmm.Sigma(muIS);
        
            outputVector{6} = muS;
            outputVector{7} = sigmaS;
            
            gmmCP = gmm.ComponentProportion(muIS);
            g = gmdistribution(muS,sigmaS);
            G = @(x)gmmCP*(g.pdf(x));
        
            [muR,muIR] = max(gmm.mu);
            sigmaR = gmm.Sigma(muIR);
            gmmCPR = gmm.ComponentProportion(muIR);
            g2 = gmdistribution(muR,sigmaR);
            G2 = @(x)gmmCPR*(g2.pdf(x));
        
            if PDFsave
                plot(finelDosages',f(finelDosages'),'-b','linewidth',2,'DisplayName',bugs{j});
                hold on
                plot(finelDosages',G2(finelDosages'),'-r','linewidth',2,'DisplayName','R');
                plot(finelDosages',G(finelDosages'),'-k','linewidth',2,'DisplayName','S');
                if ~isnan(XECOFF)
                    plot(XECOFF,0,'or','linewidth',1,'markersize',16,'DisplayName','X-ECOFF');
                end
            end
    
            tc = TCnumbers(j);
            if ~isnan(tc)
                outputVector{16} = tc;
                outputVector{17} = CIL(j);        
                outputVector{18} = CIR(j);
                if PDFsave
                    plot([CIL(j) CIR(j)],[0 0],'-r','linewidth',5,'DisplayName','T-ECOFF\pmCI');        
                    plot(tc,0,'.r','markersize',40,'HandleVisibility','off');
                end
            end
            if PDFsave
                legend boxoff
                legend('location','northwest')
                axis tight
                xlabel('log2 MIC dose')
                ylabel('likelihood')
            end
            if PDFsave && j == 1
                title(thisDrugSp)
            end
    
            if PDFsave
                subplot(1,3,3);
            end
    
            RLmatrix = exp(-abs((AIC - AIC')));
            Lmatrix = log10(RLmatrix - eye(Ngmmodels));
            matrixMax = max(max(Lmatrix));
            outputVector{15} = matrixMax;
        
            if PDFsave
                imagesc(RLmatrix)
        
                set(gca,'Xtick',1:Ngmmodels)
                set(gca,'Ytick',1:Ngmmodels)
                for k1 = 1:Ngmmodels
                    for k2 = (k1+1):Ngmmodels
                        text(k1-0.25,k2,num2str(log10(RLmatrix(k1,k2)),3),'Color',[1 1 1])
                    end
                end
                colorbar
            end
        
            if gmm.Converged
                outputVector{4} = 1;
                switch I
                    case 1
                        outputVector{25} = gmm.ComponentProportion(1);
                    case 2
                        outputVector{25} = gmm.ComponentProportion(1);
                        outputVector{27} = gmm.ComponentProportion(2);
                    case 3
                        outputVector{25} = gmm.ComponentProportion(1);
                        outputVector{26} = gmm.ComponentProportion(2);
                        outputVector{27} = gmm.ComponentProportion(3);
                end
                if PDFsave
                    title(['best is ',num2str(I),' clusters'])
                end
            else
                outputVector{4} = 0;
                if PDFsave
                    title(['best is ',num2str(I),' clusters, GMM not conv.'])
                end
            end
        
            if I > 1
                %subplot(1,4,4);
                %Pst = @(x)posterior(gmm,x);
                %plot(finelDosages,Pst(finelDosages'),'-')
                %x = solveForWTBoundary(lDosages,Pst);
                %hold on
                %plot(max(x),0,'.r','markersize',24)
                %axis tight
                %xlim([-10 10])
                %xlabel('log2 MIC dose')
                %ylabel('posterior probability')

                outputVector{10} = muR;
                outputVector{11} = sigmaR;
        
                SRboundary = fsolve(@(z)(log2(G(z))-log2(G2(z))),(muS+muR)/2,foptions);
                if PDFsave
                    subplot(1,3,2);
                    axis tight
                    yL = ylim;
                end
                if I > 2
                    [Muu,MI] = mink(gmm.mu,2);
                    sigmaInt = gmm.Sigma(MI(2));
                    muInt = gmm.mu(MI(2));
                    gmmCPInt = gmm.ComponentProportion(MI(2));
                    gInt = gmdistribution(muInt,sigmaInt);
                    GInt = @(x)gmmCPInt*(gInt.pdf(x));        
                    SIboundary = fsolve(@(x)(log10(G(x))-log10(GInt(x))),muInt,foptions);
                    IRboundary = fsolve(@(x)(log10(G2(x))-log10(GInt(x))),(muInt+muR)/2,foptions);                

                    outputVector{8} = muInt;
                    outputVector{9} = sigmaInt;
                    outputVector{12} = SIboundary;
                    outputVector{14} = IRboundary;

                    if PDFsave
                        plot([SIboundary SIboundary],[yL(1) yL(2)],'--b','DisplayName','SI boundary')    
                        plot([IRboundary IRboundary],[yL(1) yL(2)],'--r','DisplayName','IR boundary')
                    end
                end
                if PDFsave
                    plot([SRboundary SRboundary],[yL(1) yL(2)],'--k','DisplayName','SR boundary')
                end
                outputVector{13} = SRboundary;
            end
            
            bugName = split(bugs{j},' ');
            A = bugName{1};
            B = bugName{2};
        
            A = [A,'----'];
            B = [B,'----'];
            
            if length(bugName) > 2
                C = bugName{3};
                C = [C,'----'];
                bugName = [A(1:4),'-',B(1:4),C(1:2)];
            else
                bugName = [A(1:4),'-',B(1:4)];
            end
        
            if PDFsave
                export_fig(['./figures/',thisDrugSp,'_',bugName,'.pdf'])
            end
    
            [outputCell{j+2,:}] = deal(outputVector{:});
            [entireSIRData{entireCount,:}] = deal(outputVector{:});
            entirePDFs{entireCount} = f;
            entireCount = entireCount + 1;
    
            if ~isnan(CLSIS) || ~isnan(tc)
                [consolidatedSIRData{consolidatedCount,:}] = deal(outputVector{:});
                consolidatedCount = consolidatedCount + 1;
            end
        end
        
        outputTable = cell2table(outputCell(2:end,:),'VariableNames',outputCell(1,:));
        writetable(outputTable,['./csv/',thisDrugSp,'.csv'],'Delimiter',',');
    
    end
end

disp(['We have saved ',num2str(sum(~isnan(cell2mat(entireSIRData(3:end,20))))),' CLSI breakpoints'])

entireTable = cell2table(entireSIRData(2:end,:),'VariableNames',entireSIRData(1,:));
writetable(entireTable,'./csv/_entireDataset.csv','Delimiter',',');

consolidatedTable = cell2table(consolidatedSIRData(2:end,:),'VariableNames',consolidatedSIRData(1,:));
writetable(consolidatedTable,'./csv/_consolidatedDataset.csv','Delimiter',',');

save('./mat/extractedEUCAST-SIR-Data.mat','consolidatedSIRData','entireSIRData','entirePDFs')

%%

warning('on')
makeData_EUCASTToMat;


