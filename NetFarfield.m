netFarfield = outMU(:,:,1,2)- outTaper(:,:,1,2);
figure;
imagesc(abs(netFarfield));
colorbar;