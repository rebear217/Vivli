function compareANNEecDecisions(eANNdecisionDatatable,cANNdecisionDatatable)
clc
close all

[en,em] = size(eANNdecisionDatatable);
[cn,cm] = size(cANNdecisionDatatable);

realECDiffS = NaN(cn,1);
anneECDiffS = NaN(cn,1);
realECDiffR = NaN(cn,1);
anneECDiffR = NaN(cn,1);

DiffeS = NaN(cn,1);
DiffeR = NaN(cn,1);
DiffcS = NaN(cn,1);
DiffcR = NaN(cn,1);

EUCASTbpS = eANNdecisionDatatable.EUCASTbpS;
EUCASTbpR = eANNdecisionDatatable.EUCASTbpR;
CLSIbpS = cANNdecisionDatatable.CLSIbpS;
CLSIbpR = cANNdecisionDatatable.CLSIbpR;
eANNbpS = eANNdecisionDatatable.ANNbpS;
eANNbpR = eANNdecisionDatatable.ANNbpR;
cANNbpS = cANNdecisionDatatable.ANNbpS;
cANNbpR = cANNdecisionDatatable.ANNbpR;

EUCASTbpS(EUCASTbpS < -9) = NaN;
EUCASTbpR(EUCASTbpS < -9) = NaN;

for j = 1:cn
    fjBugs = find(string(cANNdecisionDatatable.bugs{j}) == string(eANNdecisionDatatable.bugs));
    fjDrugs = find(string(cANNdecisionDatatable.drugs{j}) == string(eANNdecisionDatatable.drugs));
    fj = intersect(fjDrugs,fjBugs);
    if ~isempty(fj)
        realECDiffS(fj) = EUCASTbpS(fj) - CLSIbpS(fj);
        anneECDiffS(fj) = eANNbpS(fj) - cANNbpS(fj);        
        realECDiffR(fj) = EUCASTbpR(fj) - CLSIbpR(fj);
        anneECDiffR(fj) = eANNbpR(fj) - cANNbpR(fj);
        
        DiffeS(fj) = EUCASTbpS(fj) - eANNbpS(fj);
        DiffeR(fj) = EUCASTbpR(fj) - eANNbpR(fj);
        DiffcS(fj) = CLSIbpS(fj) - cANNbpS(fj);
        DiffcR(fj) = CLSIbpR(fj) - cANNbpR(fj);
    end
end

L = 6;

figure(1)
set(1,'pos',[569         393        1251         495])

subplot(1,2,1)
plot([-L L],[-L L],'-k','linewidth',2)
hold on
plot(realECDiffS,anneECDiffS,'.b','markersize',22)
axis tight
grid on
xlabel('real S breakpoint \Delta: EUCAST - CLSI')
ylabel('ANNE S breakpoint \Delta: EUCAST - CLSI')

subplot(1,2,2)
plot([-L L],[-L L],'-k','linewidth',2)
hold on
plot(realECDiffR,anneECDiffR,'.r','markersize',22)
axis tight
grid on
xlabel('real R breakpoint \Delta: EUCAST - CLSI')
ylabel('ANNE R breakpoint \Delta: EUCAST - CLSI')

figure(2)
set(2,'pos',[569         393        650         495])

X = [-L L];

plot([-L L],[-L L],'-k','linewidth',2,'DisplayName','x=y')
hold on
plot(DiffeS,DiffeR,'.k','markersize',22,'DisplayName','(\DeltaS,\DeltaR) data')
axis tight
grid on
xlabel('S breakpoint \Delta: EUCAST - ANNE_e')
ylabel('R prediction \Delta: EUCAST - ANNE_e')
legend('location','northwest')

J = isnan(DiffeS) | isnan(DiffeR);
be = deming(DiffeS(~J),DiffeR(~J));
rhoE = corr(DiffeS(~J),DiffeR(~J));

fe = @(x)(be(1) + be(2)*x);
plot(X,fe(X),'--k','DisplayName',['Deming regression (\rho \approx ',num2str(rhoE,2),')']);

figure(3)
set(3,'pos',[569         393        650         495])

plot([-L L],[-L L],'-k','linewidth',2,'DisplayName','x=y')
hold on
plot(DiffcS,DiffcR,'.k','markersize',22,'DisplayName','(\DeltaS,\DeltaR) data')
axis tight
grid on
xlabel('S breakpoint \Delta: CLSI - ANNE_c')
ylabel('R prediction \Delta: CLSI - ANNE_c')
legend('location','northwest')

J = isnan(DiffcS) | isnan(DiffcR);
bc = deming(DiffcS(~J),DiffcR(~J));
rhoC = corr(DiffcS(~J),DiffcR(~J));

fc = @(x)(bc(1) + bc(2)*x);
plot(X,fc(X),'--k','DisplayName',['Deming regression (\rho \approx ',num2str(rhoC,2),')']);

figure(4)
set(4,'pos',[569         393        600         495])

BinWidth = 0.075;
X = -L:BinWidth:L;
X = X(:);

[Hedata,gmmE] = P(DiffeS,DiffeR);
[Hcdata,gmmC] = P(DiffcS,DiffcR);

histogram(Hedata,'Normalization','probability','BinWidth', BinWidth,'DisplayName','SR breakpoint \Delta: EUCAST - ANNE_e')
hold on
histogram(Hcdata,'Normalization','probability','BinWidth', BinWidth,'DisplayName','SR breakpoint \Delta: CLSI - ANNE_c')
axis tight
xlim([-2 2])

legend('location','northwest');
legend('boxoff')
xlabel('breakpoint \Delta (MIC)')

Ye = gmmE.pdf(X);
Ye = Ye/sum(Ye);
Yc = gmmC.pdf(X);
Yc = Yc/sum(Yc);

plot(X,Ye,'-b','DisplayName','EUCAST GMM');
plot(X,Yc,'-r','DisplayName','CLSI GMM');

col = {[0.85,0.45,0.45],[0.45 0.45 0.85],[0.45,0.85,0.45]};

figure(2)
labelled = 0;
for j = 1:length(DiffeR)
    J = gmmE.cluster(Hedata(j));
    if ~isnan(J)
        if labelled
            plot(DiffeS(j),DiffeR(j),'o','color',col{J},'linewidth',1,'markersize',9,'HandleVisibility','off');
        else
            plot(DiffeS(j),DiffeR(j),'o','color',col{J},'linewidth',1,'markersize',9,'DisplayName','cluster markers');
            labelled = 1;
        end
    end
end

figure(3)
labelled = 0;
for j = 1:length(DiffcR)
    J = gmmC.cluster(Hcdata(j));
    if ~isnan(J)
        if labelled
            plot(DiffcS(j),DiffcR(j),'o','color',col{J},'linewidth',1,'markersize',9,'HandleVisibility','off');
        else
            plot(DiffcS(j),DiffcR(j),'o','color',col{J},'linewidth',1,'markersize',9,'DisplayName','cluster markers');
            labelled = 1;
        end
    end
end

end

function [z,gmm] = P(x,y)
    %X = (X,U)U + (X,V)V
    %|X|^2 = |(X,U)|^2 + |(X,V)|^2
    x = x(:);
    y = y(:);
    X = [x , y];
    O = ones(size(X))/sqrt(2);
    O(:,2) = -O(:,2);
    U = sum(X.*O,2);
    z = U;

    MI = 200;
    options = statset('MaxIter',MI);
    gmfit{1} = fitgmdist(z,1,'Options',options,'Replicates',20,'RegularizationValue',0.01);
    gmfit{2} = fitgmdist(z,2,'Options',options,'Replicates',20,'RegularizationValue',0.01);
    gmfit{3} = fitgmdist(z,3,'Options',options,'Replicates',20,'RegularizationValue',0.01);

    AICcrit = [gmfit{1}.AIC/gmfit{1}.Converged, gmfit{2}.AIC/gmfit{2}.Converged , gmfit{3}.AIC/gmfit{3}.Converged]
    [~,jm] = min(AICcrit);
    gmm = gmfit{jm};
end
