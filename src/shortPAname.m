function str = shortPAname(drug,bug)

    drugSplit = split(drug,' ');
    bugSplit = split(bug,' ');
    len = 4;

    str = bugSplit{1};
    N = length(str);
    str = str(1:min(len,N));
    for j = 2:length(bugSplit)
        s = bugSplit{j};
        N = length(s);
        s = s(1:min(len,N));
        str = [str,s,' '];
    end
    str = str(1:end-1);
    str = [str,'-'];

    for j = 1:length(drugSplit)
        s = drugSplit{j};
        N = length(s);
        s = s(1:min(len,N));
        str = [str,s,' '];
    end

end