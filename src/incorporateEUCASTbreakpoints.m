function SIRData = incorporateEUCASTbreakpoints(SIRData,bothBreakpointSets)

    [n,~] = size(SIRData);
    for j = 3:n
        drug = SIRData{j,1};
        bug = SIRData{j,2};

        nameChange = 0;
        
        bugSplit = split(bug,' ');
        Head = lower(bugSplit{1});
        Tail = bugSplit{2:end};

        isenterobacter = 0;
        list = Enterobacterales();
        for jj = 1:length(list)
            isenterobacter = isenterobacter || strcmpi(Head,list{jj});
        end

        if isenterobacter 
            %Head = 'Enterobacterales ';
            %bug = join(Head,Tail);
            bug = 'Enterobacterales';
            nameChange = 1;
        end

        if strcmpi(Head,'Acinetobacter')
            bug = 'Acinetobacter';
            nameChange = 1;
        end

        if strcmpi(Head,'Bacteroides')
            bug = 'Bacteroides';
            nameChange = 1;
        end

        if strcmpi(Head,'Prevotella')
            bug = 'Prevotella';
            nameChange = 1;
        end

        if strcmpi(Head,'Pseudomonas')
            bug = 'Pseudomonas';
            nameChange = 1;
        end

        if strcmpi(Head,'Staphylococcus')
            %if this is Staph but not Staph aur then change to Staph...
            if ~strcmpi(Tail(1:3),'aur')
                bug = 'Staphylococcus';
                nameChange = 1;
            end
        end
        
        %match SIR data name to either CLSI or EUCAST:
        D = upper(drug);
        B = upper(bug);

        Fd = find( (D == bothBreakpointSets.Antibiotic_eucast) | (D == bothBreakpointSets.Antibiotic_clsi));
        Fb = find( (B == bothBreakpointSets.Bacteria_eucast) | (B == bothBreakpointSets.Bacteria_clsi));
        F = intersect(Fd,Fb);

        %now check for a good match or errors and correct them if possible:
        %A good match:
        if ~isempty(F) && (length(F) == 1)
            SIRData{j,23} = bothBreakpointSets.S_eucast(F);
            SIRData{j,24} = bothBreakpointSets.R_eucast(F);
        end
        %No match:
        if isempty(F)
            disp(['No match found for ',SIRData{j,1},' & ',SIRData{j,2}])
            if nameChange
                disp(['It did not help to match against ',D,' & ',B])            
            end
        end

        %several matches, one of which could be wrong:
        if length(F) > 1
            if (length(unique(bothBreakpointSets.S_eucast(F))) == 1) && (length(unique(bothBreakpointSets.R_eucast(F))) == 1)
                %there a actually is only 1 breakpoint here, so we can use it:
                SIRData{j,23} = bothBreakpointSets.S_eucast(F(1));
                SIRData{j,24} = bothBreakpointSets.R_eucast(F(1));
            else
                %if there are several breakpoints in the matched list, we can attemtpt to resolve them by matching again within this matched list...
                disp(['An incorrect matching issue was found with ',D,' & ',B,' : ...'])
                disp(['Antibiotic eucast',',','Antibiotic clsi',',','Bacterium eucast',',','Bacterium clsi',',','S bp eucast',',','R bp eucast'])
    
                for ii = 1:length(F)
                    Aeucast = bothBreakpointSets.Antibiotic_eucast(F(ii));
                    Aclsi = bothBreakpointSets.Antibiotic_clsi(F(ii));
    
                    Beucast = bothBreakpointSets.Bacteria_eucast(F(ii));
                    Bclsi = bothBreakpointSets.Bacteria_clsi(F(ii));
    
                    Seucast = bothBreakpointSets.S_eucast(F(ii));
                    Reucast = bothBreakpointSets.R_eucast(F(ii));
    
                    disp([Aeucast,',',Aclsi,',',Beucast,',',Bclsi,',',num2str(Seucast),',',num2str(Reucast)])
                end
                %let's find just the EUCAST matches of drug and bug in the matched already matched list
                %(indexed by F):
                Fde = find(D == bothBreakpointSets.Antibiotic_eucast(F));
                Fbe = find(B == bothBreakpointSets.Bacteria_eucast(F));
                %this is drug and bug matches, let's make sure we can find a unique one:
                Fe = intersect(Fde,Fbe);
                attemptResolveBP_S = unique(bothBreakpointSets.S_eucast(F(Fe)));
                attemptResolveBP_R = unique(bothBreakpointSets.R_eucast(F(Fe)));

                %if the S and R EUCAST breakpoints actually just have 1
                %breakpoint in them, then we're ok and we can correct the multi-matches now by honing in on the right (and unique) breakpoint:
                if (length(attemptResolveBP_S) == 1) && (length(attemptResolveBP_R) == 1)
                    SIRData{j,23} = attemptResolveBP_S;
                    SIRData{j,24} = attemptResolveBP_R;
                    disp('The error was resolved because there is a unique correct match within the (broken) set of all possible matches.')
                    disp(' ')
                else
                    %at this point, we still have many matches for some
                    %reason:
                    disp('The error could not be resolved.')
                    disp(' ')
                end

            end

        end

        %if nameChange
        %   disp(['No match for ',drug,' & ',oldbug,' : also tried matching as ',bug])
        %else
        %   disp(['No match for ',drug,' & ',oldbug])
        %end

    end

end