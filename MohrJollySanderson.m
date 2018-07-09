%   MohrJollySanderson.m - script to plot Mohr circles, with J&S opening
%   fractures
%   
%   equations from:
%       Jolly & Sanderson, 1997 Journal of Structural Geology
%
%   David Healy 
%   August 2016 
%   d.healy@abdn.ac.uk

close all ; 

disp(' ') ; 
disp(['*** Started MohrJollySanderson.m at ', datestr(now), '...']) ; 

disp(' ') ; 
disp('Fractures...') ; 
%fnFractures = input('Filename of fracture pole data: ') ; 
fnFractures = 'sarahspoles.txt' ; 
fidFractures = fopen(fnFractures, 'r') ; 
[polesFractures, nFractures] = fscanf(fidFractures, '%g %g', [2, inf]) ; 
fclose(fidFractures) ; 
nFractures = nFractures / 2 ; 
polesFractures = polesFractures' ; 
polesFracturesRad = polesFractures * pi / 180 ; 

%   cartesian/Andersonian stresses in MPa
disp(' ') ; 
disp('Stresses...') ; 
%sigmaV = input('Enter vertical stress, MPa: ') ; 
%sigmaHmax = input('Enter maximum horizontal stress, MPa: ') ; 
%sigmaHmin = input('Enter minimum horizontal stress, MPa: ') ; 
sigmaV = 75 ; 
sigmaHmax = 50 ; 
sigmaHmin = 30 ; 
sorted_sigma = sort([ sigmaV sigmaHmax sigmaHmin ]) ; 
pf = 60 ; 

%   failure envelope 
mu = 0.6 ; 
C0 = 20 ; 
sigmaN = 100 ; 
tau = C0 + mu * sigmaN ; 

%   principal stresses
sigma1 = sorted_sigma(3) ; 
sigma2 = sorted_sigma(2) ; 
sigma3 = sorted_sigma(1) ; 
sigmad = sigma1 - sigma3 ; 

Phi = ( sigma2 - sigma3 ) / ( sigma1 - sigma3 ) ; 
Rprime = ( pf - sigma3 ) / ( sigma1 - sigma3 ) ; 

disp(' ') ; 
disp(['Stress ratio \Phi: ', num2str(Phi)]) ; 
disp(['Stress ratio R'': ', num2str(Rprime)]) ; 

%   calculate normal and shear stresses on supplied fractures     
for f = 1:nFractures 

    phi_rad = polesFracturesRad(f, 1) ; 
    theta_rad = polesFracturesRad(f, 2) ;
        
    %   direction cosines
    l = sin(theta_rad) * cos(phi_rad) ; 
    m = cos(theta_rad) * cos(phi_rad) ; 
    n = -sin(phi_rad) ; 

    %   normal and tau stress on plane with this pole 
    fsigmaN(f) = ( sigmaHmin * l^2 ) + ( sigmaHmax * m^2 ) + ( sigmaV * n^2 ) ; 
    ftau(f) = sqrt( ( ( sigmaHmin - sigmaHmax )^2 * l^2 * m^2 ) + ... 
                   ( ( sigmaHmax - sigmaV )^2 * m^2 * n^2 ) + ...
                   ( ( sigmaV - sigmaHmin )^2 * n^2 * l^2 ) ) ;
        
end ; 

%   plot the Mohr diagram, with fractures  
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 6 6 ]) ; 

hold on ; 
plotMohr(sigma1, sigma2, sigma3, 1, '-k') ; 
plot([0, sigmaN], [C0, tau], '-r', 'LineWidth', 1) ; 
fill([0, pf, pf, 0], [0, 0, sigma1, sigma1], 'b', 'FaceAlpha', 0.2) ; 
for f = 1:nFractures
    plot(fsigmaN(f), ftau(f), '.k', 'MarkerSize', 10) ; 
end ; 
hold off ; 
xlim([0 sigma1*1.05]) ; 
ylim([0 sigmad*1.1]) ;
xlabel('Effective normal stress, MPa') ; 
ylabel('Shear stress, MPa') ; 
title({'Mohr circles, after Jolly & Sanderson';''}) ; 

print -r300 -djpeg 'MohrJollySanderson.jpeg' ; 

disp(' ') ; 
disp(['*** ...finished MohrJollySanderson.m at ', datestr(now), '.']) ; 
disp(' ') ; 
