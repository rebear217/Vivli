function [bothBreakpointSetsOut,SIRDataOut] = validateSIRdata(bothBreakpointSetsIn,SIRDataIn)

data = extractSIRdataColumns(SIRDataIn);
[n,m] = size(bothBreakpointSetsIn);
[~,N] = size(data.drugs);

count = 0;
errorCount = 0;
updateCount = 0;

for j = 1:N
    drug = upper(data.drugs{j});
    bug = upper(data.bugs{j});
    PApair = [bug,' ',drug];
    F = find(upper(PApair) == bothBreakpointSetsIn.FullName_clsi);
    if ~isempty(F)
        count = count + 1;

        A = bothBreakpointSetsIn.S_clsi(F);
        B = bothBreakpointSetsIn.R_clsi(F);
        C = bothBreakpointSetsIn.S_eucast(F);
        D = bothBreakpointSetsIn.R_eucast(F);
        
        a = data.CLSIsticS(j);
        b = data.CLSIsticR(j);
        c = data.eucastBPS(j);
        d = data.eucastBPR(j);
    
        V = [A B C D];
        v = [a b c d];
        Dif = norm(v-V,1);
        if Dif > 0
            errorCount = errorCount + 1;
        end
        if isnan(Dif)
            disp(['Found potential variable differences for ',PApair])
            makeNote = 1;
            if isnan(a) && ~isnan(A)
                SIRDataIn{j+2,19} = bothBreakpointSetsIn.S_clsi(F);
                updateCount = updateCount + 1;
                if makeNote
                    disp('Updated variable')
                    makeNote = 0;
                end
            end
            if isnan(b) && ~isnan(B)
                SIRDataIn{j+2,20} = bothBreakpointSetsIn.R_clsi(F);
                updateCount = updateCount + 1;
                if makeNote
                    disp('Updated variable for ')
                    makeNote = 0;
                end
            end
            if isnan(c) && ~isnan(C)
                SIRDataIn{j+2,23} = bothBreakpointSetsIn.S_eucast(F);
                updateCount = updateCount + 1;
                if makeNote
                    disp('Updated variable')
                    makeNote = 0;
                end
            end
            if isnan(d) && ~isnan(D)
                SIRDataIn{j+2,24} = bothBreakpointSetsIn.R_eucast(F);
                updateCount = updateCount + 1;
                if makeNote
                    disp('Updated variable')
                    makeNote = 0;
                end
            end
            if isnan(A) && ~isnan(a)
                bothBreakpointSetsIn.S_clsi(F) = SIRDataIn{j+2,19};                
                updateCount = updateCount + 1;
                if makeNote
                    disp('Updated variable')
                    makeNote = 0;
                end
            end
            if isnan(B) && ~isnan(b)
                bothBreakpointSetsIn.R_clsi(F) = SIRDataIn{j+2,20};
                updateCount = updateCount + 1;
                if makeNote
                    disp('Updated variable for ')
                    makeNote = 0;
                end
            end
            if isnan(C) && ~isnan(c)
                bothBreakpointSetsIn.S_eucast(F) = SIRDataIn{j+2,23};
                updateCount = updateCount + 1;
                if makeNote
                    disp('Updated variable')
                    makeNote = 0;
                end
            end
            if isnan(D) && ~isnan(d)
                bothBreakpointSetsIn.R_eucast(F) = SIRDataIn{j+2,24};
                updateCount = updateCount + 1;
                if makeNote
                    disp('Updated variable')
                    makeNote = 0;
                end
            end
            if makeNote
                disp('No point updating variables')
            end
        end
        
    end
    
end

disp(['Checked ',num2str(count),' data entries']);
disp(['Found ',num2str(errorCount),' errors where CLSI & EUCAST S/R differed between variables']);
disp(['Updated ',num2str(updateCount),' data entries due to NaNs']);

SIRDataOut = SIRDataIn;
bothBreakpointSetsOut = bothBreakpointSetsIn;

end