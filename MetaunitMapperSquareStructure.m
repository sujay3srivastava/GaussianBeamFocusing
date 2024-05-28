clear;

% This code takes the square surface dimensions and find the center point of each 
% meta unit//

wgx= 100e-6;  % waveguide rectangle length in x
wgy= 100e-6; % waveguide rectangle length in y
ay = 0.5e-6; % metaunit property ~lambda_eff/2
ax = 0.45e-6; % ax not equal to ay but close
ys = 0;
ye = wgy;

rad = 0.15e-6; % radius of metaunits
y = 0:2*ay:wgy;
Ny = length(y);
Coordinates = {};
for ny=1:Ny-1  
    Y0 = y(ny)+ay; % iterating through different points at distance 1 metaunit away in y direction
    if (Y0+ay+(2*rad))<ye

        X0 = -0.5*wgx + 2*ax; % centerline
        Coordinates{end+1} = [X0,Y0];
        nx = 1;
        while wgx/2 > (X0+2*ax+ ax) %Different metaunits creation where metaunits are populated in a square
%area
            X0= X0+2*ax;
            Coordinates{end+1} = [X0,Y0];
            nx = nx + 1;

        end
    end
end
     
matrixCoordinates = vertcat(Coordinates{:}); % separates (x,y) to  adjacent elements in same row

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds the target beam, performs the mode correction///
% Example usage:
% Set parameters
% Function to calculate 2D Gaussian intensity
%load("phasecolormap.mat");
%load("deltacolormap.mat");
A = 1;                      % Amplitude
x0 = 0;                     % x-coordinate of the center
y0 = (wgy)/2;             % y-coordinate of the center
z = 1e-4;                   % Distance from focal point
lambda= 1.55e-6;            % Wavelength
k = 2*pi/lambda;            % wave number
W0 =  24e-6;                % beam width at center
z0 = W0^2/lambda*pi;        % Rayleigh Range
Wz = W0*sqrt(1+(z/z0)^2);
Rz= z*(1+ (z0/z)^2);
zeta = atan(z/z0);

%matrix Coordinates contains points where the Intensity/Phase need to be
%calculated
%Finding required Phase and Intensity for gaussian beam focus
xpoints = linspace(-wgx/2, wgx/2,2 *1e6 * (wgx) / 0.90);
ypoints = linspace(0,wgy,wgy*1e6);
[P,Q] =meshgrid(xpoints,ypoints);
intensity_test = gaussian2D(P,Q, A, x0, y0, W0,Wz);
phase_test = gaussianphase2D(P,Q, A, x0, y0, k,Rz,z,zeta);
% The above is for making figure on taper profile

f = figure;
width= 200;
height=200;
f.Position(3:4) = [width height];
figure(f);
imagesc(intensity_test);
colorbar;
xlabel('X');
ylabel('Y');
zlabel('Intensity Target');
title('2D Gaussian Intensity Function');


% These function are used for finding gaussian profile at metaunit points
intensity = gaussian2D(matrixCoordinates(:, 1),matrixCoordinates(:, 2), A, x0, y0, W0,Wz);
phase = gaussianphase2D(matrixCoordinates(:, 1),matrixCoordinates(:, 2), A, x0, y0, k,Rz,z,zeta);

figure;
imagesc(phase_test);
colorbar;
xlabel('X');
ylabel('Y');
zlabel('Phase Target');
%Getting mode correction data from taper simulation
load('F:\MATLABCodesLMW\GaussianBeam_DeltaCode\GaussianBeam-SquareStructure\SquareTaperFields\square_field_100um.mat');
figure;
Ez  = Ez(:,:,1,1);
imagesc(x_E, y_E, transpose(abs(Ez)));
colorbar; 
Ez_interp = interp2(x_E, y_E, transpose(Ez), matrixCoordinates(:, 1),matrixCoordinates(:, 2));
% figure; 
% colorbar;
% imagesc(matrixCoordinates(:, 1),matrixCoordinates(:, 2), abs(Ez_interp));
% title('E-field Amplitude on the Taper');


%Performing Mode Correction on Gaussian Intensity and Phase
maxEz = max(max(abs(Ez_interp)));
Corrected_amplitude = intensity./abs(Ez_interp/maxEz);

% Compensating Phase  (Changing to count individual point phase)
% nX = zeros(numel(matrixCoordinates(:, 1)),1);
% centralPhase = angle(interp2(x_E, y_E, transpose(Ez), nX,matrixCoordinates(:, 2)));
% compensationPhase = angle(Ez_interp) - centralPhase;
% Corrected_phase = phaseCor(phase - compensationPhase);

compensationPhase = angle(Ez_interp);
Corrected_phase = phaseCor(phase - compensationPhase);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code used to find delta values for each meta unit, the data will be
%passed to lumerical script to create the meta-units


numRows = size(matrixCoordinates, 1);  %total points
Delta1 = zeros(numRows, 1);
Delta2 = zeros(numRows, 1);

for i = 1:numRows
   
     Delta1(i,1)= del1cal(Corrected_amplitude(i), Corrected_phase(i));
     Delta2(i,1) = del2cal(Corrected_amplitude(i), Corrected_phase(i));
    
end
% scaling the values to cap the max delta to 200 nm, hence we find maximum
% delta value and scale accordingly

scale1 = max(abs(Delta1));  
scale2 = max(abs(Delta2));
scale = max(scale1,scale2);
Delta1 = (Delta1/scale)*200e-9;
Delta2 = (Delta2/scale)*200e-9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DELTA-SCALING%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FUNCTIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function intensity = gaussian2D(x, y, A, x0, y0, W0,Wz)
    % 2D Gaussian function for amplitude
    intensity = A * W0 / Wz* exp(-1 *((x - x0).^2 + (y - y0).^2 )/Wz^2);
end

function phase = gaussianphase2D(x, y, ~, x0, y0, k, Rz,z,zeta)
    % 2D Gaussian function for phase
    phase = k*z+ k*(0.5*((x - x0).^2  + (y - y0).^2)/Rz) - zeta;
    phase = phase - 2*pi*floor(phase/2/pi +0.5);           % making range [-pi,pi]
end

function phi = phaseCor(phi)
    [Nx, Ny] = size(phi);
    for nx = 1:Nx
        for ny = 1:Ny
            if abs(phi(nx,ny)) > pi
                phi(nx,ny) = phi(nx,ny) - sign(phi(nx,ny))*2*pi;
            end
        end
    end
end

function delta1 = del1cal(Corrected_amplitude, Corrected_phase)
    delta1 = abs(Corrected_amplitude).*cos(Corrected_phase);
end
function delta2 = del2cal(Corrected_amplitude, Corrected_phase)
    delta2 = abs(Corrected_amplitude).*sin(Corrected_phase);
end

