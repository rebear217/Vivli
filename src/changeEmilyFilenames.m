clc
d = dir('EUCAST_vs_ATLAS/histograms/*.pdf');

for j = 1:length(d)
    oldname = d(j).name;
    newname = strrep(oldname,'_','-');
    newname = strrep(newname,' ','-');
    %pause
    if ~strcmp(oldname,newname)
        disp([num2str(j),' : ',oldname,' : ',newname])
        movefile(['./EUCAST_vs_ATLAS/histograms/',oldname],...
            ['./EUCAST_vs_ATLAS/histograms/',newname])
    end
end

%%

disp(' ')
disp(' ')
disp(' ')

wildcards = {'*levo*.pdf','*mino*.pdf','*ceftaz*.pdf','*dapto*.pdf'};
disp('------------')
for k = 1:length(wildcards)
    d = dir(['EUCAST_vs_ATLAS/histograms/',wildcards{k}]);
    for j = 1:length(d)
        disp(d(j).name)
    end
    disp('------------')
end

