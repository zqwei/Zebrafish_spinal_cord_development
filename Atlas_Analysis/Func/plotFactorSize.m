function plotFactorSize(factorSize)
nNeurons = size(factorSize, 1);
timePoints = (1:size(factorSize, 2))*240;
nCol = 8;
nRow = ceil(nNeurons/nCol);
figure('units', 'pixels', 'outerposition', [0 0 , 1000, 1000]);
for i = 1:nNeurons
    subplot(nRow, nCol, i);
    hold on
    scatter(timePoints(factorSize(i, :)>2)/3600/4, factorSize(i, factorSize(i, :)>2), 5, 'b');
    scatter(timePoints(factorSize(i, :)==1)/3600/4, factorSize(i, factorSize(i, :)==1), 5, 'r');
    scatter(timePoints(factorSize(i, :)==2)/3600/4, factorSize(i, factorSize(i, :)==2), 5, 'g');
    xlim([0 5])
    hold off
end

for i = 1:nNeurons
    subplot(nRow, nCol, i);
    hold on
    scatter(timePoints(factorSize(i, :)>2)/3600/4, factorSize(i, factorSize(i, :)>2), 5, 'b');
    scatter(timePoints(factorSize(i, :)==1)/3600/4, factorSize(i, factorSize(i, :)==1), 5, 'r');
    scatter(timePoints(factorSize(i, :)==2)/3600/4, factorSize(i, factorSize(i, :)==2), 5, 'g');
    ylim([0 10]);
    xlim([0 3]);
    hold off
end