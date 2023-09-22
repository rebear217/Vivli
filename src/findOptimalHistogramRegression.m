function [params,fit,linearFit] = findOptimalHistogramRegression(YvectorData,XmatrixData,initialguess)

    %Suppose histograms are laid out with different drugs vertically and different
    %dosages horizonally, so say histogramsC/E = 55 (down) x 19 (across) and
    %then breakpointDiff = 55 column vector. We then want to find the weight
    %vector W that maximises the correlation between histogramsC - histogramsE
    %and breakpointDiff in some norm. The purpose of this code is to
    %construct that norm.
    
    [~,m] = size(XmatrixData);
    
    d = YvectorData(:);
    D = XmatrixData;

    if nargin < 3
        initialguess = ones(m+1,1);
        initialguess(1) = 0;
        initialguess = initialguess/sum(initialguess);
        initialguess = 0*initialguess;
    end

    %The idea is that d = model(p,histogramsC - histogramsE) for some p:
    
    opts = statset('Display','iter','TolFun',1e-12,'MaxIter',1000,'MaxFunEvals',1000,'Jacobian','on');
    %opts = statset();

    fitPlus = fitnlm(D,d,@(p,D)model(p,D,+1),initialguess,'Options',opts);
    fitMinus = fitnlm(D,d,@(p,D)model(p,D,-1),initialguess,'Options',opts);

    %There are 2 possible fits depending on whether the correlation is + or
    %-ve and this copes with the constraints that weights >= 0 and add up
    %to unity.

    %Thus the model d = a + b(DW) where sum(W) = 1, W >= 0 (a and b real numbers)
    %is re-written as a pair of unconstrained fits:
    %d = a + sD(w.^2) where s = +1 or -1. Then b = s*sum(w.^2) and 
    %W = w.^2 / sum(w.^2).

    %So now choose the best fit:
    %RL = exp(-abs(fitMinus.ModelCriterion.AIC - fitPlus.ModelCriterion.AIC));

    if fitPlus.Rsquared.Ordinary > fitMinus.Rsquared.Ordinary
        fit = fitPlus;
        S = +1;
    else
        fit = fitMinus;
        S = -1;
    end

    a = fit.Coefficients.Estimate(1);
    w2 = abs(fit.Coefficients.Estimate(2:end));
    b = S*sum(w2);
    W = w2/sum(w2);

    %Error bounds:
    %da and d(w^2) are reported from fitnlm, dS = 0, abs(s) = 1:
    %da = fit.Coefficients.SE(1);
    %dw2 = fit.Coefficients.SE(2:end);

    %db = d(sum(w^2)) = sum(d(w^2)): want +ve value for this
    %db = abs(sum(dw2));

    %|dW| = d(w^2) / sum(w^2) + w^2*db/(sum(w^2))^2: using + not - sign
    %from the deriative to get upper bound on SE:
    %dW = (dw2 / b) - w2*db/b^2;

    %params.RL = RL;
    params.a = a;
    params.b = b;
    params.w = W;
    %params.da = da;
    %params.db = db;
    %params.dw = dw2;

    weightedModel = a + b*(XmatrixData*W);

    linearFit.YvectorData = YvectorData;
    linearFit.XweightedFit = weightedModel;
    linearFit.correlation = corr(weightedModel,YvectorData);
    
end

function R = model(p,D,mysign)

    %Parameter vector p holds the y-intercept (a) and weights (w):
    a = p(1);
    w = p(2:end);%typically 19-dimensional

    %We want the constraint that w >=0 so we square the input weight to make it
    %positive:
    
    w = w.^2;
    
    %We now need to weight each dose by w by forming w(1)D(:,1) + w(2)D(:,2) + ... + w(55)D(:,55) = Dw
    %which is of the shape 55x19 * 19x1 = 55x1. Then we add a to make the
    %linear model for the breakpoint difference based on the histogram
    %ensemble difference:
    
    R = a + mysign*D*w;
    
    %We can now try and fit this R against d data above

end
