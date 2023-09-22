function mergeValidatePriorNewANNdecisions(ANNdecisionData,newHistosANNdecisionData,CEflag)

    close all

    if nargin < 3
        CEflag = 1;
    end    
    
    [n,~] = size(newHistosANNdecisionData.table);
    Fs = [];
    Js = [];

    for j = 1:n
        drug = newHistosANNdecisionData.table.drugs{j};
        bug = newHistosANNdecisionData.table.bugs{j};

        Fd = find(strcmpi(drug,ANNdecisionData.table.drugs));
        Fb = find(strcmpi(bug,ANNdecisionData.table.bugs));
        F = intersect(Fd,Fb);
        if ~isempty(F)
            if length(F) > 1
                error(['Multiple matches issued occurred with ',drug,' and ',bug])
            else
                if CEflag
                    newHistosANNdecisionData.table.CLSIbpS(j) = ANNdecisionData.table.CLSIbpS(F);
                    newHistosANNdecisionData.table.CLSIbpR(j) = ANNdecisionData.table.CLSIbpR(F);
                else
                    newHistosANNdecisionData.table.EUCASTbpS(j) = ANNdecisionData.table.EUCASTbpS(F);
                    newHistosANNdecisionData.table.EUCASTbpR(j) = ANNdecisionData.table.EUCASTbpR(F);
                end
                Fs = [Fs F];
                Js = [Js j];
            end
        end

    end

    disp('The figures should look the same, with wider IQR in Figure 2 to reflect the linear ANN "correction" step ...')
    ANNanalyseOneHistogram('Ceftazidime','Klebsiella pneumoniae',ANNdecisionData,CEflag,1)
    ANNanalyseOneHistogram('Ceftazidime','Klebsiella pneumoniae',newHistosANNdecisionData,CEflag,2)

    dotColour = [1 1 1]/2;
    compareANNwithBP(newHistosANNdecisionData,CEflag,3,[1 1 1]/2);
    figure(3)
    title('corrected ANN decisions')
    figure(4)
    title('corrected ANN decisions')
    compareANNwithBP(ANNdecisionData,CEflag,5);
    figure(5)
    title('prior ANN decisions')
    figure(6)
    title('prior ANN decisions')

    X = [-7 7];

    figure(7)
    [rho,pval] = corr(ANNdecisionData.table.ANNbpS(Fs),newHistosANNdecisionData.table.ANNbpS(Js),'type','Pearson');
    plot(ANNdecisionData.table.ANNbpS(Fs),newHistosANNdecisionData.table.ANNbpS(Js),'.','color',dotColour,'DisplayName',['prior/post ANN decisions check (r \approx ',num2str(rho,3),')'],'MarkerSize',22)
    hold on
    plot(X,X,'-k','DisplayName','x=y')
    axis tight
    legend()

    figure(8)
    [rho,pval] = corr(ANNdecisionData.table.ANNbpR(Fs),newHistosANNdecisionData.table.ANNbpR(Js),'type','Pearson');
    plot(ANNdecisionData.table.ANNbpR(Fs),newHistosANNdecisionData.table.ANNbpR(Js),'.','color',dotColour,'DisplayName',['prior/post ANN decisions check (r \approx ',num2str(rho,3),')'],'MarkerSize',22)
    hold on
    plot(X,X,'-k','DisplayName','x=y')
    axis tight
    legend()
    
    if CEflag
        CnewHistosANNdecisionData = newHistosANNdecisionData;
        save('./mat/c_allnewHistosANNdecisions.mat','CnewHistosANNdecisionData')
    else
        EnewHistosANNdecisionData = newHistosANNdecisionData;
        save('./mat/e_allnewHistosANNdecisions.mat','EnewHistosANNdecisionData')
    end

end