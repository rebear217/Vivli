function [list,output,check] = Enterobacterales(input)

    if nargin < 1
        input = 'XXXXX';
    end

    s = split(input,' ');
    s = s{1};
    input = lower(s);

    list = {'enterobacterales',...
        'enterobacter',...
        'escherichia',...
        'salmonella',...
        'klebsiella',...
        'citrobacter',...
        'morganella',...
        'serratia',...
        'proteus',...
        'providencia',...
        'yersinia'};
    
    output = input;
    check = 0;

    for l = 1:length(list)
        L = list{l};
        if strcmpi(L,input)
            output = 'enterobacterales';
            check = 1;
        end
    end

end