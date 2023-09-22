
disp(['total CLSI breakpoints ', ...
    num2str(sum(~isnan(bothBreakpointSets.R_clsi)))]);

disp(['total EUCAST human-curated (stringent) text and breakpoint matches ', ...
    num2str(sum(~isnan(bothBreakpointSets.R_eucast)))])

disp(['total CLSI PAs that yield a EUCAST PA match (relaxed) by sceening the CLSI PA list ', ...
    num2str(sum(~isnan(cell2mat(consolidatedSIRData(3:end,19)))))])

disp(['total CLSI PAs with a text EUCAST PA match (relaxed) by screening the CLSI PA list ' ...
    'that also have a valid EUCAST breakpoint ',num2str(sum(~isnan(cell2mat(consolidatedSIRData(3:end,23)))))])

disp(['total relaxed CLSIs-EUCAST matches where both breakpoints exist ',...
    num2str(sum(~isnan(cell2mat(consolidatedSIRData(3:end,23)) .* cell2mat(consolidatedSIRData(3:end,19)))))])

%%
