addpath('../Data_Analysis/')
for nFile = [6 11 24]
    FAEV_v0_1(nFile);
    Leader_v0(nFile);
    Leader_v1_2(nFile);
    Lineage_v1(nFile);
end