function numbers = TECOFFStoNumbers(TECOFFS)

    numbers = NaN(size(TECOFFS));
    
    for n = 1:length(TECOFFS)
        tc = string(TECOFFS(n));
        if strcmp(tc(1),'-') && (length(tc) > 1)
            tc = tc(2:end);
        end
        if ~ismember(tc,{'ID','-'})
            tcnum = NaN;
            try
                tcnum = str2num(tc);
            catch
                if strcmp(tc(1),'(')
                    if strcmp(tc(2),'(')
                        tc = tc(3:end-2);
                        tcnum = str2num(tc);
                    else
                        tc = tc(2:end-1);
                        tcnum = str2num(tc);
                    end
                end
            end
            if tcnum < 0
                tcnum = -tcnum;
            end
            numbers(n) = tcnum;
        end
    end
end