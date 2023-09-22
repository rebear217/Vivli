function allCLSIEUCASTMatches = findAllCLSIEUCASTMatches(SIRData,EUCASTBreakpoints)
%This file is no longer used: REB 19 Feb 2023

    %SIRData = consolidatedSIRData;
    %SIRData = entireSIRData;
    
    [n,~] = size(SIRData);
    allCLSIEUCASTMatches = cell(n-1,9);
    allCLSIEUCASTMatches{1,1} = 'CLSI Abx';
    allCLSIEUCASTMatches{1,2} = 'EUCAST Abx';
    allCLSIEUCASTMatches{1,3} = 'CLSI bacterium';
    allCLSIEUCASTMatches{1,4} = 'EUCAST bacterium';
    allCLSIEUCASTMatches{1,5} = 'another possible match';
    allCLSIEUCASTMatches{1,6} = 'CLSI BP S';
    allCLSIEUCASTMatches{1,7} = 'CLSI BP R';
    allCLSIEUCASTMatches{1,8} = 'EUCAST BP S';
    allCLSIEUCASTMatches{1,9} = 'EUCAST BP R';

    for j = 2:n-1
        if mod(j,50) == 0 
            disp(j/n)
        end
        thisCLSIDrug = SIRData{j+1,1};
        thisCLSIBug = SIRData{j+1,2};
        [matches,~,~,foundBugs,foundDrugs] = matchOneCLSItoEUCAST(thisCLSIDrug,thisCLSIBug,EUCASTBreakpoints);
        allCLSIEUCASTMatches{j,1} = thisCLSIDrug;
        allCLSIEUCASTMatches{j,3} = thisCLSIBug;
        allCLSIEUCASTMatches{j,6} = SIRData{j+1,19};
        allCLSIEUCASTMatches{j,7} = SIRData{j+1,20};
        if ~isempty(matches)
            allCLSIEUCASTMatches{j,2} = foundDrugs{1};
            allCLSIEUCASTMatches{j,4} = foundBugs{1};
            allCLSIEUCASTMatches{j,8} = log2(EUCASTBreakpoints{matches(1),3});
            allCLSIEUCASTMatches{j,9} = log2(EUCASTBreakpoints{matches(1),4});

            if length(matches) > 1
                allCLSIEUCASTMatches{j,5} = [foundDrugs{2},' & ',foundBugs{2}];
            end        
        end
    end

end