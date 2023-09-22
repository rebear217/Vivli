function gmmFitsHistogramPlot(entireSIRData)

    X = cell2mat(entireSIRData(3:end,5));
    histogram(X)
    xlim([0.7 1]);
    xlabel('goodness of fit (R^2)')
    ylabel('# observations')
    
end