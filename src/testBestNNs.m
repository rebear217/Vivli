function testBestNNs(bestNets,CEflag,lDosages,testtype)

if nargin < 4
    testtype = 1;
end

close all

%P = globalParameterValues;
popSize = 1000;
%bins = lDosages;
N = 9;

s = linspace(0,1,N);
popSize1 = popSize*ones(1,N);
popSize2 = popSize*ones(1,N);
Lsigma = ones(1,N);
Rsigma = ones(1,N);
switch testtype
    case 1
        Rmean = linspace(-2,4,N);
        Lmean = -3*ones(1,N);
        popSize2 = round(popSize2/2);
        Lsigma = Lsigma/2;
        Rsigma = Rsigma/2;
    case 2
        Rmean = 3*ones(1,N);
        Lmean = -3*ones(1,N);
        popSize2 = round(popSize2.*s);
        popSize1 = round(popSize1.*(1-s));
    case 3
        Rmean = linspace(5,5,N);
        Lmean = linspace(-1,2,N);
        popSize2 = round(popSize2/2);
        Lsigma = Lsigma/2;
        Rsigma = Rsigma/2;
end

BPLs = zeros(N,1);
BPRs = zeros(N,1);
dIQRL = zeros(N,2);
dIQRR = zeros(N,2);
histoSurface = zeros(length(lDosages),N);

figure(1)
set(1,'pos',[57           1        2494        1336])

for j = 1:N
    %figure(j)
    %set(j,'pos',[228   321   816   526]);
    subplot(3,3,j)

    synthdata = syntheticData([popSize1(j) popSize2(j)],[Lmean(j) Rmean(j)],[Lsigma(j) Rsigma(j)]);
    %[~,~,~,~,ECOFF] = findECOFF(synthdata,bins,P.ECOFFpValue);

    synthdata = synthdata/sum(synthdata);
    histoSurface(:,j) = synthdata;

    ANNbreakpointDecisions = ANNdecision(synthdata',bestNets,[],CEflag);

    plot(lDosages,synthdata,'.-k','markersize',28,'DisplayName','MIC histogram');
    M = 1.4*max(synthdata);

    hold on
    
    bpL = ANNbreakpointDecisions.bpLefts;
    bpR = ANNbreakpointDecisions.bpRights;
    BPLs(j) = ANNbreakpointDecisions.bpLeftFinal;
    BPRs(j) = ANNbreakpointDecisions.bpRightFinal;
    dIQRL(j,:) = ANNbreakpointDecisions.bfLeftIQR;
    dIQRR(j,:) = ANNbreakpointDecisions.bfRightIQR;
        
    Nse = randn(size(bpL))*0.01;

    bins = min(lDosages):0.5:max(lDosages);

    histogram(bpL,bins,'FaceColor','b','Normalization','probability')
    histogram(bpR,bins,'FaceColor','r','Normalization','probability')
    
    plot(bpL,M*(0.9*ones(size(bpL)) + Nse),'.b','markersize',12,'DisplayName','bpL ANNE')
    plot([BPLs(j) BPLs(j)],[0 1]*M,'-b','DisplayName','bpL decision')
    
    plot(bpR,M*(0.95*ones(size(bpR)) + Nse),'.r','markersize',12,'DisplayName','bpR ANNE')
    plot([BPRs(j) BPRs(j)],[0 1]*M,'-r','DisplayName','bpR decision')

    plot([dIQRL(j,1) dIQRL(j,2)],M*[0.8 0.8],'-b','DisplayName','bpL IQR','linewidth',0.5)
    plot([dIQRR(j,1) dIQRR(j,2)],M*[0.85 0.85],'-r','DisplayName','bpR IQR','linewidth',0.5)

    %plot(ECOFF,0,'ob','DisplayName','alg ECOFF','markersize',14,'linewidth',1)
    
    axis tight
    xlim([lDosages(1) lDosages(end)])
    set(gca,'Xtick',lDosages)
    set(gca,'XtickLabel',round(lDosages,2))
    xlabel('MIC (log2 \mug/mL)')
    ylabel('frequency (synthetic data)')
    if j == 1
        legend('location','east')
        title('ANNE breakpoint decisions: synthetic test case 1')
    else
        title(['synthetic test case ',num2str(j)])
    end

    %histogram(bpL,'Normalization','probability','BinEdges', lDosages);
    drawnow

end

fg = figure();
[X,Y] = meshgrid(1:N,lDosages);
surf(X,Y,histoSurface)
m = max(histoSurface(:));
colormap(1-bone)
hold on
for j = 1:N
    plot3(j,BPLs(j),m,'.b','markersize',34);
    plot3([j j],[BPLs(j) BPLs(j)],[0 m],'-b','LineWidth',0.5);
    plot3(j,BPRs(j),m,'.r','markersize',34);
    plot3([j j],[BPRs(j) BPRs(j)],[0 m],'-r','LineWidth',0.5);
    plot3([j j],[dIQRL(j,1) dIQRL(j,2)],[m m],'-b','LineWidth',2);
    plot3([j j],[dIQRR(j,1) dIQRR(j,2)],[m m],'-r','LineWidth',2);
end

view(2)
axis tight
set(gca,'Ytick',lDosages)
set(gca,'Yticklabel',round(lDosages,2))
xlabel('synthetic test case label')
ylabel('MIC (log2 \mug/mL)')
zlabel('frequency (synthetic)')
colorbar
