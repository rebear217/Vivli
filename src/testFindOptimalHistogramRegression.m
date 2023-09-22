
close all
clc

%%

m = 19;
n = 55;

%a random vector to represent breakpoint CLSI-EUCAST differences:
breakpointDiff = randn(n,1);
%2 random matrices to represent CLSI and EUCAST histograms
histogramsC = rand(n,m);
histogramsE = rand(n,m);
%make sure they are histograms by normalising
sC = sum(histogramsC,2);
RC = repmat(sC,1,19);
hC = histogramsC./RC;

sE = sum(histogramsE,2);
RE = repmat(sE,1,19);
hE = histogramsE./RE;

%determine the difference between histograms:
histogramsDiff = hE - hC;

%compute the optimal weight to histogram differences that makes them model breakpoint
%differences as a linear model that weights each dose:
[ps,fit,linear] = findOptimalHistogramRegression(breakpointDiff,histogramsDiff);

%Plot the results:
figure(1)
plot(1:m,ps.w,'-k');
title('weight vector')
xlabel('index j - i.e. dose level')
ylabel('w_j')
axis tight

%Now plot the linear model that makes use of the weights found above:

figure(2)
plot(breakpointDiff,linear.XweightedFit,'.k','markersize',24)
hold on
%plot(breakpointDiff,linear.feval(breakpointDiff),'-k');
title(['\rho\approx ',num2str(linear.correlation,2)])
xlabel('breakpoint difference')
ylabel('weighted-norm histogram difference')
axis tight

%%
