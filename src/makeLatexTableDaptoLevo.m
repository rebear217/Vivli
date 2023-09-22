function makeLatexTableDaptoLevo(bothBreakpointSets)

clc

FL = find(bothBreakpointSets.Antibiotic_clsi == 'LEVOFLOXACIN');
FD = find(bothBreakpointSets.Antibiotic_clsi == 'DAPTOMYCIN');

for j = 1:length(bothBreakpointSets.Antibiotic_clsi)
    display = 0;
    drug = bothBreakpointSets.Antibiotic_clsi(j);
    bug = bothBreakpointSets.Bacteria_clsi(j);
    levoSclsi = '-';
    levoRclsi = '-';
    daptoSclsi = '-';
    daptoRclsi = '-';
    levoSeucast = '-';
    levoReucast = '-';
    daptoSeucast = '-';
    daptoReucast = '-';
    if ismember(j,FL)
        display = 1;
        levoSclsi = myProcess(bothBreakpointSets.S_clsi(j));
        levoRclsi = myProcess(bothBreakpointSets.R_clsi(j));        
        levoSeucast = myProcess(bothBreakpointSets.S_eucast(j));
        levoReucast = myProcess(bothBreakpointSets.R_eucast(j));        
    end
    if display
        x = 1;
    end
    if ismember(j,FD)
        display = 1;
        daptoSclsi = myProcess(bothBreakpointSets.S_clsi(j));
        daptoRclsi = myProcess(bothBreakpointSets.R_clsi(j));        
        daptoSeucast = myProcess(bothBreakpointSets.S_eucast(j));
        daptoReucast = myProcess(bothBreakpointSets.R_eucast(j));        
    end
    if display
        S = lower(char(bug));
        S(1) = upper(S(1));
        myStr = strcat("{\it ",S,"} & ",levoSclsi,'/',levoRclsi," & ",levoSeucast,'/',levoReucast,...
            " & ",daptoSclsi,'/',daptoRclsi," & ",daptoSeucast,'/',daptoReucast," \\");
        disp(myStr);
    end
end

end

function Sout = myProcess(S)
    if isempty(S) || isnan(S)
        Sout = '-';
    else
        Sout = string(S);
    end
end
