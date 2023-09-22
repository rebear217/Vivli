function [matches,foundBugMatch,foundDrugMatch,foundBugs,foundDrugs] = matchOneCLSItoEUCAST(CLSIdrug,CLSIbug,EUCASTbreakpoints)
%This file is no longer used: REB 19 Feb 2023

    synonymsCLSI = {{'Escherichia','Salmonella',...
        'Klebsiella','Citrobacter','Enterobacter'}};
    synonymsEUCAST = {'Enterobacterales'};
    
    CLSIdrug = replace(CLSIdrug,'_',' ');
    CLSIbug = replace(CLSIbug,',',' ');

    BH = split(CLSIbug,' ');
    if length(BH) > 1
        bughead = BH{1};
        bugtail = BH{2};
        bugHandT = 1;
    else
        bughead = BH{1};
        bugtail = '';
        bugHandT = 0;
    end
    
    DH = split(CLSIdrug,' ');
    if length(DH) > 1
        drughead = DH{1};
        drugtail = DH{2};
        drugHandT = 1;
    else
        drughead = DH{1};
        drugtail = '';
        drugHandT = 0;
    end
    
    n = length(synonymsCLSI);
    for j = 1:n
        [memflag,i] = ismember(bughead,synonymsCLSI{j});
        if memflag
            bughead = synonymsEUCAST{j};
            bugtail = '';
            bugHandT = 0;
            break
        end
    end
    
    [n,~] = size(EUCASTbreakpoints);
    foundBugMatch = zeros(n,1);
    foundDrugMatch = zeros(n,1);
    
    for j = 1:n
        thisBug = EUCASTbreakpoints{j,1};
        thisDrug = EUCASTbreakpoints{j,2};
    
        thisBugSplit = split(thisBug,' ');
        if ismember(bughead,thisBugSplit)
            if ~bugHandT
                foundBugMatch(j) = 1;
            end
            if bugHandT && ismember(bugtail,thisBugSplit)
                foundBugMatch(j) = 2;
            end
        end
        thisDrug = replace(thisDrug,'-',' ');
        thisDrug = replace(thisDrug,'_',' ');        
        thisDrug = replace(thisDrug,'.','');
        thisDrug = replace(thisDrug,',','');
        thisDrug = replace(thisDrug,')','');
        thisDrug = replace(thisDrug,'(','');
        thisDrug = replace(thisDrug,'1','');
        thisDrug = replace(thisDrug,'2','');
        thisDrug = replace(thisDrug,'3','');
        thisDrug = replace(thisDrug,'4','');
        thisDrug = replace(thisDrug,'5','');
        thisDrug = replace(thisDrug,'6','');
        
        thisDrugSplit = split(thisDrug,' ');
        if ismember(drughead,thisDrugSplit)
            if ~drugHandT
                foundDrugMatch(j) = 1;
            end
            J = join(thisDrugSplit(2:end));
            J = J{1};
            J = replace(J,' ','');
            if drugHandT
                Jcheck = all(ismember(drugtail,J));
                if Jcheck
                    foundDrugMatch(j) = 2;
                    if length(drugtail) < length(J)
                        foundDrugMatch(j) = 3;
                    end
                end
            else
                if ~isempty(J)
                    foundDrugMatch(j) = 3;
                end
            end
        end
    
    end

    matches = find(foundDrugMatch&foundBugMatch);
    foundBugs = EUCASTbreakpoints(matches,1);
    foundDrugs = EUCASTbreakpoints(matches,2);
    
    if isempty(matches) && bugHandT
        [matches,foundBugMatch,foundDrugMatch,foundBugs,foundDrugs] = matchOneCLSItoEUCAST(CLSIdrug,bughead,EUCASTbreakpoints);
    end

end