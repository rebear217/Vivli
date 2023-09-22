function y = signedlog10(x)
    y = sign(x).*log10(abs(x));
end