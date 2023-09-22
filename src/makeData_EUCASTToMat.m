clear all
close all

%%

%ignoreList = {'Ceftazidime-clavulanate.csv'};
ignoreList = [];

rootDir = './EUCASThistograms/';
d = dir([rootDir,'*.csv']);

analysisFile = d(1).name;

lDosagesString = {'all EUCAST histograms','dose (log2 ug/mL)','-8.9658','-7.9658','-6.9658','-5.9658','-5.0589','-4.0589',...
    '-3','-2','-1','0','1','2','3','4','5','6','7','8','9'};
lDosages = [-8.9658,-7.9658,-6.9658,-5.9658,-5.0589,-4.0589,-3,-2,-1,0,1,2,3,4,5,6,7,8,9];

totalFileNo = length(d);
allEUCASTHistograms = cell(3,length(lDosagesString));
count = 1;

%%

for fileJ = 1:totalFileNo

    analysisFile = d(fileJ).name;
    disp(analysisFile)

    loadError = 0;
    if ismember(analysisFile,ignoreList)
        loadError = 1;
    end
    try
        MICdata = importMICs([rootDir,analysisFile]);
        bugs = importBugNames([rootDir,analysisFile]);
        %[TECOFF, ConfidenceInterval] = importECOFFs([rootDir,analysisFile]);
    catch
        disp('This has a non-standard CSV structure.')
        loadError = 1;
    end
    if isempty(bugs{1})
        loadError = 1;
        disp('This has a non-standard CSV structure (it loaded but with junk).')
    end

    if ~loadError

        thisDrugSp = getDrugFromFilename(analysisFile);

        [n,m] = size(MICdata);
        for j = 1:n
            allEUCASTHistograms{count,1} = thisDrugSp;
            allEUCASTHistograms{count,2} = bugs{j};
            outputVector = num2cell(MICdata(j,:));
            [allEUCASTHistograms{count,3:m+2}] = deal(outputVector{:});
            count = count + 1;
        end

    end

end

%%

histoTable = cell2table(allEUCASTHistograms,'VariableNames',lDosagesString);
writetable(histoTable,'./csv/_all-EUCAST-Histogtrams.csv','Delimiter',',');

save('./mat/all-EUCAST-Histograms.mat','allEUCASTHistograms','lDosagesString','lDosages')

