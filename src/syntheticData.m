function [data,bins] = syntheticData(N,mu,sigma)
    bins = -9:9;
    X = NaN(max(N),length(mu));
    for j = 1:length(mu)
        x = sigma(j)*randn(N(j),1)+mu(j);
        X(1:length(x),j) = x;
    end
    Y = X(:);
    Y = Y(~isnan(Y));
    Y = Y(Y>=min(bins));
    Y = Y(Y<=max(bins));
    data = histcounts(Y,[bins bins(end)+1]);
end