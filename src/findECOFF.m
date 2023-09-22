function [bestfit,bestbin,R2,wholeRMS,ECOFF,ECOFFCorr,allFits] = findECOFF(histogram,bins,ECOFFp)

    histogram = histogram/sum(histogram);

    jmin = find(histogram > 0,3,'first');
    N = length(bins);
    fitXdata = bins(jmin(1):jmin(end));
    fitYdata = histogram(jmin(1):jmin(end));
    model = @(p,x)gaussian(x,p(1),p(2));
    muguess = mean(fitXdata);
    sigmaguess = 1;
    fit = fitnlm(fitXdata,fitYdata,model,[muguess sigmaguess]);
    wholeRMS = NaN(1,N);
    R2 = NaN(1,N);
    R2(jmin(end)) = fit.Rsquared.Ordinary;
    bestfit = fit;
    bestbin = jmin(end);
    bestR2 = R2(jmin(end));
    
    allFits = cell(1,N);
    allFits{jmin(end)} = fit;

    for b = (jmin(end)+1):N
        fitXdata = bins(jmin(1):b);
        fitYdata = histogram(jmin(1):b);
        muguess = fit.Coefficients.Estimate(1);
        sigmaguess = fit.Coefficients.Estimate(2);    
        fit = fitnlm(fitXdata,fitYdata,model,[muguess sigmaguess]);
        allFits{b} = fit;

        R2(b) = fit.Rsquared.Ordinary;
        if R2(b) > bestR2
            bestfit = fit;
            bestbin = b;
            bestR2 = R2(b);
        end
        wholeRMS(b) = norm(fit.feval(bins)-histogram,2);
    end
    allFits = {allFits{(jmin(end)):N}};

    mu = bestfit.Coefficients.Estimate(1);
    sigma = bestfit.Coefficients.Estimate(2);
    ECOFF = icdf('Normal',ECOFFp,mu,sigma);

    F = find(bins < ECOFF);
    H = histogram(F);
    H = H / norm(H,2);
    M = bestfit.feval(bins(F));
    M = M / norm(M,2);
    ECOFFCorr = dot(H,M);
    
    if ECOFF > max(bins) || ECOFF < min(bins) || abs(imag(ECOFF)) > 1e-15
        ECOFF = NaN;
        R2 = NaN;
        bestbin = NaN;
        wholeRMS = NaN;
        ECOFFCorr = NaN;
    end
end
