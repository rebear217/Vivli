function CLSIenterobacteralesCount(bothBreakpointSets,rawCLSIbreakpoints)
    
    clc
    P = globalParameterValues();
    lDosages = P.lDosages;

    for fig = 1:2
        if fig == 1
            drugs = bothBreakpointSets.Antibiotic_clsi;
            bugs = bothBreakpointSets.Bacteria_clsi;
            breakPs_S = bothBreakpointSets.S_clsi;
            breakPs_R = bothBreakpointSets.R_clsi;
        else
            drugs = rawCLSIbreakpoints.Antibiotic;
            bugs = rawCLSIbreakpoints.Species;
            breakPs_S = rawCLSIbreakpoints.S;
            breakPs_R = rawCLSIbreakpoints.R;
        end
    
        udrgS = unique(drugs);
        ubugS = unique(bugs);
        
        count = 0;
        Ecount = 0;
    
        DataS = {};
        DataR = {};
        Drugs = {};
    
        for i = 1:length(udrgS)
    
            thisDrug = udrgS(i);
            disp('----------------------------------------------')
            disp(['-- ',thisDrug,' --'])
            disp('----------------------------------------------')
    
            DataS{i} = [];
            DataR{i} = [];
            Drugs{i} = string(thisDrug);
    
            fS = find(drugs == thisDrug);
            for j = 1:length(ubugS)
                
                thisBug = ubugS(j);
                [~,~,Echeck] = Enterobacterales(string(thisBug));
                fR = find(bugs == thisBug);
                F = intersect(fR,fS);
                if ~isempty(F)
                    for f = 1:length(F)
                        if Echeck
                            disp([bugs(F(f)), ...
                                num2str(breakPs_S(F(f))), ...
                                num2str(breakPs_R(F(f)))])
                            Ecount = Ecount + 1;
                            DataS{i} = [DataS{i} breakPs_S(F(f))];
                            DataR{i} = [DataR{i} breakPs_R(F(f))];
                        end
                        count = count + 1;
                    end
                end
            end
        end
    
        disp(['Checked ',num2str(count),' entries in the CLSI list'])
        disp(['Found ',num2str(Ecount),' enterobacterales in the CLSI list'])
    
        figure(fig)
        set(fig,'pos',[157         816        1466         412])
    
        for i = 1:length(Drugs)
            bpS = DataS{i};
            bpR = DataR{i};
            if ~isempty(bpS)
                if length(unique(bpS)) == 1
                    marker = '*k';
                    ms = 12;
                    p = plot(i*ones(size(bpS)) + 0*(rand(size(bpS))-0.5),bpS,marker,'markersize',ms);
                    text(i+0.1,bpS(1)+0.25,num2str(length(bpS)));
                else
                    marker = '.k';            
                    ms = 30;
                    P = plot(i*ones(size(bpS)) + 0*(rand(size(bpS))-0.5),bpS,marker,'markersize',ms);
                end
            end
            hold on
            if ~isempty(bpR)
                if length(unique(bpR)) == 1
                    marker = '*r';
                    ms = 12;
                    q = plot(i*ones(size(bpR)) + 0*(rand(size(bpR))-0.5),bpR,marker,'markersize',ms);
                    text(i+0.1,bpR(1)+0.25,num2str(length(bpR)),'color','r');
                else
                    marker = 'sr';            
                    ms = 12;
                    Q = plot(i*ones(size(bpS)) + 0*(rand(size(bpS))-0.5),bpR,marker,'markersize',ms);
                end
            end
    
            %{
            if length(unique(bpS)) > 1
                errorbar(i+0.3,mean(bpS),sqrt(var(bpS)),'ok','markersize',8);
            end
            if length(unique(bpR)) > 1
                errorbar(i+0.6,mean(bpR),sqrt(var(bpR)),'or','markersize',8);
            end
            %}
        end
    
        set(gca,'Xtick',1:length(Drugs))    
        set(gca,'XtickLabels',Drugs)
        set(gca,'Ytick',round(lDosages,2))    
        set(gca,'YtickLabels',round(lDosages,2))
        ylabel('MIC (log2 \mug/mL)')
    
        legend([p,q,P,Q],{'unique CLSI S breakpoints','unique CLSI R breakpoints','non-unique CLSI S breakpoints','non-unique CLSI R breakpoints'})
        legend('location','southeast')
        title([num2str(Ecount),'/',num2str(count),' CLSI breakpoints are for Enterobacterales strains'])
    
    end
end