clear all
close all
clc
warning('off')

%%

Params = globalParameterValues();

PDFsave = 0;
%ignoreList = {'Ceftazidime-clavulanate.csv'};
ignoreList = [];

breakpoints = importCLSIbreakpoint('./CLSI/breakpoints.csv');
ECOFFpValue = Params.ECOFFpValue;

disp(['Imported ',num2str(length(breakpoints.S)),' breakpoints'])

outputTitles = {'antibiotic','pathogen','optimal cluster number','GMM converged','R squared',...
    'S mean','S sigma','I mean','I sigma','R mean','R sigma',...
    'SI boundary','SR boundary','IR boundary','max log10 relative likelihood',...
    '(T)ECOFF','(T)ECOFF CI min','(T)ECOFF CI max','CLSI STIC S','CLSI STIC R',...
    '(X)ECOFF','(X)ECOFF correlation','EUCAST S bp','EUCAST R bp','GMMweightS','GMMweightI','GMMweightR'};
outputUnits = {'string','string','integer','Boolean','real',...
    'log2 ug/mL','(log2 ug/mL)^2','log2 ug/mL','(log2 ug/mL)^2','log2 ug/mL','(log2 ug/mL)^2',...
    'log2 ug/mL','log2 ug/mL','log2 ug/mL','real',...
    'log2 ug/mL','log2 ug/mL','log2 ug/mL','log2 ug/mL','log2 ug/mL',...
    'log2 ug/mL','real','log2 ug/mL','log2 ug/mL','real','real','real'};

%%

rootDir = './EUCASThistograms/';
d = dir([rootDir,'*.csv']);

analysisFile = d(1).name;

Dosages = importDosages([rootDir,analysisFile]);
lDosages = log2(Dosages);
finelDosages = min(lDosages):0.01:max(lDosages);
fineDosages = 2.^finelDosages;

totalFileNo = length(d);
entireSIRData = cell(3,length(outputTitles));
entireCount = 3;
consolidatedSIRData = cell(3,length(outputTitles));
consolidatedCount = 3;

entirePDFs = {};

nancount = 0;
okcount = 0;
matchList = {};

for fileJ = 1:totalFileNo
    
    close all

    analysisFile = d(fileJ).name;
    disp('--------------------------');
    disp(analysisFile)
    disp('--------------------------');

    loadError = 0;
    if ismember(analysisFile,ignoreList)
        loadError = 1;
    end
    try
        MICdata = importMICs([rootDir,analysisFile]);
        bugs = importBugNames([rootDir,analysisFile]);
        [TECOFF, ConfidenceInterval] = importECOFFs([rootDir,analysisFile]);
    catch
        disp('This has a non-standard CSV structure.')
        loadError = 1;
    end
    if isempty(bugs{1})
        loadError = 1;
        disp('This has a non-standard CSV structure (it loaded but with junk).')
    end

    if ~loadError
        TCnumbers = log2(TECOFFStoNumbers(TECOFF));
        [CIL,CIR] = CItoNumbers(ConfidenceInterval);
        CIL = log2(CIL);
        CIR = log2(CIR);
        
        %%

        thisDrugSp = getDrugFromFilename(analysisFile);
        disp(['    Check : ',analysisFile,' provides data for this drug ... ',thisDrugSp])

        %%
        
        [n,m] = size(MICdata);
    
        %%
        for j = 1:n
            [CLSIS,CLISR,spA,spB,abA,abB] = getBreakpointNEW(bugs{j},thisDrugSp,breakpoints);
            if isnan(CLSIS)
                nancount = nancount + 1;
            else
                [memcheck , loc] = ismember([spA,spB,abA,abB],matchList);
                if memcheck
                    disp(['Prior match found : ',bugs{j},',',thisDrugSp,',',spA,',',spB,',',abA,',',abB])
                    %disp(['Prior match found : ',matchList{loc}]);
                end
                okcount = okcount + 1;
                matchList{okcount} = [spA,spB,abA,abB];
            end
        end
        disp(['Matches now up to ',num2str(okcount)]);
    end
end

matchList = matchList';

disp('--------------------------------------------------')
disp(['Imported ',num2str(length(breakpoints.S)),' breakpoints'])
disp(['Found ',num2str(okcount),' breakpoints'])
disp(['Discarded ',num2str(nancount),' attempted matches'])
disp('--------------------------------------------------')
