function y = gaussian(x,mu,sigma)
    x = x(:);
    y = exp(-(x-mu).^2/2/sigma.^2)/sigma/sqrt(2*pi);
end