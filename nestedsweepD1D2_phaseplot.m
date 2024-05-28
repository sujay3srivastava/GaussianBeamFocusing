% load("AmplitudeSweepX45degreemesh4.mat");
load("phasecolormap.mat");
%bias =  by(1,1,6,6);
%cal_by = by - bias;
%dim_values = abs(squeeze(cal_by(1, 1, :, :))); % Extracting the values from 'by'
dim_values = angle(squeeze(bx(1, 1, :, :))); %- by(1,1,6,6))); % Extracting the values from 'by'
bx_squeeze = (squeeze(bx(1, 1, :, :)))
figure();
imagesc(delta1,delta2, dim_values);
camroll(90);
c = phasecolourmap;
colormap(c);
colorbar;
xlabel("delta1");
ylabel("delta2");
title("Phase Graph")