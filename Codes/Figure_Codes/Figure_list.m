% figure 2_d
control_datasets = [3, 4, 7, 12, 10, 11, 13, 15, 16];

h = figure('Position', [0, 0, 1500, 200]);
hold on;
for nFile = control_datasets
    Figure_2(h, nFile);
end
