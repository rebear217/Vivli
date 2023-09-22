function DeOptimiseQ(weight,cbestNets,ebestNets)

    close all

    P = globalParameterValues();
    P.lDosages = round(P.lDosages,2);
    bins = min(P.lDosages):0.33:max(P.lDosages);
    
    options = optimset('Display','iter','MaxFunEvals',20);
    
    CEflags = [0 0 1 1];
    SorRflags = [0 1 0 1];

    for fCount = 1:4
        H0 = [0 0 0 0 0 0 5 25 50 200 50 25 5 0 0 0 0 0 0];
        H0 = H0/sum(H0);

        CEflag = CEflags(fCount);
        SorRflag = SorRflags(fCount);

        figure(fCount)
        set(fCount,'Position',[125   816   666   402]);

        if CEflag
            bestNets = cbestNets;
            bp0 = ANNdecision(H0',cbestNets,[],CEflag);
            str = 'CLSI';
        else
            bestNets = ebestNets;
            bp0 = ANNdecision(H0',ebestNets,[],CEflag);
            str = 'EUCAST';
        end

        if SorRflag
            str = [str,' S'];
        else
            str = [str,' R'];
        end

        itLimit = 80;
        it = 0;
        BPdifference = 0;
        Hguess = H0;
        while (BPdifference < 1) && (it < itLimit)
            Hguess = fminsearch(@(H)Q(H),Hguess',options)';
            Hguess = abs(Hguess);
            Hguess = Hguess/sum(Hguess);
            bp = ANNdecision(Hguess',bestNets,[],CEflag);
            if SorRflag
                BPdifference = abs(bp.bpLeftFinal - bp0.bpLeftFinal)
            else
                BPdifference = abs(bp.bpRightFinal - bp0.bpRightFinal)
            end
            it = it + 1;
        end

        plot(P.lDosages,H0,'.-b','linewidth',1,'DisplayName','MIC distribution start (h_0)','markersize',24)
        hold on
        plot(P.lDosages,Hguess,'.-r','linewidth',1,'DisplayName','MIC distribution end (h)','markersize',24)

        if SorRflag
            histogram(bp0.bpLefts,bins,'Normalization','probability','FaceColor','b','DisplayName','ANNE start')
            histogram(bp.bpLefts,bins,'Normalization','probability','FaceColor','r','DisplayName','ANNE end')

            plot([bp0.bpLeftFinal bp0.bpLeftFinal],[0 1],'--b','linewidth',1,'DisplayName',[str,' breakpoint start'])
            plot([bp.bpLeftFinal bp.bpLeftFinal],[0 1],'--r','linewidth',1,'DisplayName',[str,' breakpoint end'])
        else
            histogram(bp0.bpRights,bins,'Normalization','probability','FaceColor','b','DisplayName','ANNE start')
            histogram(bp.bpRights,bins,'Normalization','probability','FaceColor','r','DisplayName','ANNE end')

            plot([bp0.bpRightFinal bp0.bpRightFinal],[0 1],'--b','linewidth',1,'DisplayName',[str,' breakpoint start'])
            plot([bp.bpRightFinal bp.bpRightFinal],[0 1],'--r','linewidth',1,'DisplayName',[str,' breakpoint end'])
        end

        set(gca,'Xtick',P.lDosages)
        ylabel('frequency')
        xlabel('MIC (log2 \mug/mL)')
        ylim([0 0.6])
        xlim([-6 8])
        legend('location','northwest')
        legend('boxoff')
        drawnow;

    end

    function q = Q(H)

        H = abs(H);
        H = H/sum(H);
        bp = ANNdecision(H,bestNets,[],CEflag);
        if SorRflag
            q = norm(weight.*(H-H0)) / abs(bp.bpLeftFinal - bp0.bpLeftFinal);
        else
            q = norm(weight.*(H-H0)) / abs(bp.bpRightFinal - bp0.bpRightFinal);
        end
    
    end

end

