function drugname = getDrugFromFilename(analysisFile)

        thisDrug = split(analysisFile,'.');
        thisDrug = thisDrug{1};
        thisDrugSp = split(thisDrug,'-');
        if length(thisDrugSp) > 1
            thisDrugSp = [thisDrugSp{1},' ',thisDrugSp{2}];
        else
            thisDrugSp = thisDrug;
        end
        if strcmp(thisDrugSp(end),' ')
            thisDrugSp = thisDrugSp(1:end-1);
        end
  
        drugname = thisDrugSp;
        drugname = replace(drugname,'-',' ');
        drugname = replace(drugname,'_',' ');
end