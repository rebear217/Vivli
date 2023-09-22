function summariseHistograms(allEUCASTHistograms)

data = cell2mat(allEUCASTHistograms(:,3:end));
[n,m] = size(data);

for j = 1:n
    data(j,:) = data(j,:)/sum(data(j,:));
end

figure(1)
imagesc(data);
view(2);
colormap('jet')
set(gca,'Xtick',1:20)
set(gca,'Xticklabel',-10:10)
xlabel('MIC (log2 ug/mL)')
ylabel('Pathogen-antibiotic pair')
colorbar

end