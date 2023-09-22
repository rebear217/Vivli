function testCountBreakpointMatchings(consolidatedSIRData,bothBreakpointSets)

    [n,~] = size(bothBreakpointSets);
    consolidatedSIRData = incorporateEUCASTbreakpoints(consolidatedSIRData,bothBreakpointSets);
    data = extractSIRdataColumns(consolidatedSIRData);
    disp('-------')
    disp(['There are ',num2str(sum(~isnan(data.eucastBPS))),' breakpoint S entries in consolidatedSIRData.']);
    disp(['There are ',num2str(sum(~isnan(data.eucastBPR))),' breakpoint R entries in consolidatedSIRData.']);
    disp(['Out of ',num2str(length(data.eucastBPR)),' total entries.']);
    disp(['There are ',num2str(n),' valid CLSI S and R breakpoints.']);
    disp('-------')

end