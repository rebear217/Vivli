function bpHistograms(bothBreakpointSets)

    close all
    bins = -9:1:9;
    figure(1)
    set(1,'pos',[98         823        1104         407])

    subplot(1,2,1)
    histogram(bothBreakpointSets.S_clsi,'FaceColor','r','BinEdges',bins,'Normalization','probability','DisplayName','breakpoint S CLSI')
    hold on
    histogram(bothBreakpointSets.S_eucast,'FaceColor','b','BinEdges',bins,'Normalization','probability','DisplayName','breakpoint S EUCAST')
    legend('location','northwest')
    xlabel('MIC log2 \mug/mL')
    ylabel('frequency')

    subplot(1,2,2)
    histogram(bothBreakpointSets.R_clsi,'FaceColor','r','BinEdges',bins,'Normalization','probability','DisplayName','breakpoint R CLSI')
    hold on
    histogram(bothBreakpointSets.R_eucast,'FaceColor','b','BinEdges',bins,'Normalization','probability','DisplayName','breakpoint R EUCAST')
    legend('location','northwest')
    xlabel('MIC log2 \mug/mL')
    ylabel('frequency')

end