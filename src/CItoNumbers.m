function [CInumbersL,CInumbersR] = CItoNumbers(ConfidenceInterval)

    CInumbersL = NaN(size(ConfidenceInterval));
    CInumbersR = NaN(size(ConfidenceInterval));
    
    for n = 1:length(ConfidenceInterval)
        if ~isempty(ConfidenceInterval{n})
            s = split(ConfidenceInterval{n},'-');
            CInumbersL(n) = str2num(s{1});
            CInumbersR(n) = str2num(s{2});
        end
    end

end