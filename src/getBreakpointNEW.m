function [S,R,spAout,spBout,abAout,abBout] = getBreakpointNEW(bug,antibiotic,breakpoints)
%the antibiotic variable needs spaces between words in it, not '-' signs, to work properly

    antibiotic = char(antibiotic);
    bug = char(bug);
    
    antibiotic = replace(antibiotic,'-',' ');
    antibiotic = replace(antibiotic,'_',' ');
    
    spA = shorten(char(bug));
    abA = shorten(char(antibiotic));

    N = length(breakpoints.Species);
    same = false;

    j = 1;
    foundJ = NaN;
    
    S = NaN;
    R = NaN;

    while ~same && (j <= N)
        spB = shorten(char(breakpoints.Species(j)));
        abB = shorten(char(breakpoints.Antibiotic(j)));

        sameBug = strcmpi(spA,spB);
        sameAb = strcmpi(abA,abB);
        same = sameAb && sameBug;

        if same
            foundJ = j;
        end

        j = j + 1;
    end

    if ~isnan(foundJ)
        S = breakpoints.S(foundJ);
        R = breakpoints.R(foundJ);
        spAout = spA;
        spBout = spB;
        abAout = abA;
        abBout = abB;
    else
        spAout = '';
        spBout = '';
        abAout = '';
        abBout = '';
    end
end

function shortstring = shorten(longstring)
    SB = split(longstring,' ');
    n = length(SB);
    sb1 = SB{1};
    sb1 = [sb1,'******'];
    if n > 1
        sb2 = SB{2};
        sb2 = [sb2,'******'];
        if n >= 2
            shortstring = [sb1(1:6),sb2(1:6)];
        else
            %this is never reached: 3-string lables are truncted to 2
            %labels because CLSI only use 2 whereas EUCAST differentiat
            %this into 3. But this may lead to "over-matching":
            sb3 = SB{3};
            sb3 = [sb3,'******'];
            shortstring = [sb1(1:5),sb2(1:5),sb3(1:5)];
        end
    else
        shortstring = sb1(1:6);
    end

end