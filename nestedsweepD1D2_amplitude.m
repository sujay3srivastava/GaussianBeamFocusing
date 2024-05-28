%load("r150delta50mesh6.mat");
load("deltacolormap.mat");
%bias =  by(1,1,6,6);
%cal_by = by - bias;
%dim_values = abs(squeeze(cal_by(1, 1, :, :))); % Extracting the values from 'by'
dim_values = abs(squeeze(by(1, 1, :, :))); % Extracting the values from 'by'
figure();
imagesc(delta2,delta1, dim_values);
camroll(90)
c = deltacolormap;
colormap(c);
colorbar;

xlabel("delta2");
ylabel("delta1");
title("Amplitude")
