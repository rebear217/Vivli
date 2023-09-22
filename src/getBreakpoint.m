function [S,R,spAout,spBout,abAout,abBout] = getBreakpoint(bug,antibiotic,breakpoints)

%the antibiotic variable needs spaces between words in it, not '-' signs, to work properly
warning('getBreakpoint does not work properly, it gives too many matches, use getBreakpointNEW instead.')

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
    sb1 = [sb1,'-----'];
    if n > 1
        sb2 = SB{2};
        sb2 = [sb2,'---'];
        if n == 2
            shortstring = [sb1(1:5),sb2(1:3)];
        else
            sb3 = SB{3};
            sb3 = [sb3,'---'];
            shortstring = [sb1(1:5),sb2(1:3),sb3(1:3)];
        end
    else
        shortstring = sb1(1:5);
    end
    %if strcmpi(shortstring,'Amp') && (n > 1)
    %    shortstring = longstring(1:4);
    %end
end