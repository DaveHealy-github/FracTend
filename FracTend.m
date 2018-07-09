%   FracTend.m - script to plot slip and dilatation tendency  
%   
%   Equations & code from:
%       Morris et al., 1996 Geology
%       Ferrill et al., 1999 GSA Today 
%       Streit & Hillis, 2004 Energy
%       Jolly & Sanderson, 1997 Journal of Structural Geology 
%       Allmendinger et al., 2012 Structural Geology Algorithms, Cambridge
%       University Press
%
%   David Healy & Tara Stephens  
%   July 2018 
%   d.healy@abdn.ac.uk

close all ; 

version = 0.9 ; 

disp(' ') ; 
disp(['*** Started FracTend version ', num2str(version), ' at ', datestr(now), '...']) ; 
disp(' ') ; 

%   read in poles to specific fractures; tab-delimited text file, 
%   formatted as plunge then trend 
disp(' ') ; 
disp('Reading input file of poles to planes...') ; 
%fnFractures = input('Filename of fracture pole data: ') ; 
fnFractures = 'Utah_OA-Sills.txt' ;  % 'Utah_Thrusts_TS-RW. Utah_Sills_ArcRW. Utah_DefBands2RW. Utah_OA-Sills
fidFractures = fopen(fnFractures, 'r') ; 
[polesFractures, nFractures] = fscanf(fidFractures, '%g %g', [2, inf]) ; 
fclose(fidFractures) ; 
nFractures = nFractures / 2 ; 
polesFractures = polesFractures' ; 
polesFracturesRad = polesFractures * pi / 180 ;
disp(['Read ', num2str(nFractures), ' fracture poles']) ; 

%   read in stress magnitudes 
%   principal stresses in MPa
disp(' ') ; 
disp('Stresses...') ; 
sigma1 = 43 ;      
sigma2 = 38.86 ;      
sigma3 = 25 ;       
sorted_sigma = [ sigma1, sigma2, sigma3 ] ; 
sigmad = sigma1 - sigma3 ; 
disp(['Principal stresses ', num2str(sorted_sigma), ' in MPa']) ; 

% pore fluid pressure in MPa, for fracture Opening Angle
% calculations 
Pf = 37.24;

%   note the implicit convention: x//s1, y//s2, z//s3
stressTensor = [ sorted_sigma(1), 0, 0 ; ...
                 0, sorted_sigma(2), 0 ; ... 
                 0, 0, sorted_sigma(3) ] ; 

%   read in stress orientation 
%   e.g. for old case of SHmax azimuth of 135 
        % normal fault system: Trend s1=0, Plunge s1=90 - use s3-trend to change orientation of stress field with s1 as the rotation axis
        % thrust fault: Trend of s3 must be >90 from s1; Ps1=0
        % strike-slip: Trend of s3 must be 90 from Trend of s1
trendS1 = 68 ; 
plungeS1 =  3 ;    
trendS3 = 265;   
disp('Stress orientation:') ; 
disp(['sigma1 plunge/trend - ', num2str(plungeS1, '%02d'), '/', num2str(trendS1, '%03d')]) ;                     
disp(['sigma3 trend - ', num2str(trendS3, '%03d')]) ;                     
trendS1rad = trendS1 * pi / 180 ; 
plungeS1rad = plungeS1 * pi / 180 ; 
trendS3rad = trendS3 * pi / 180 ; 

%   coefficient of friction & cohesion 
disp(' ') ; 
muStatic = 0.6 ; 
C0 = 0 ; 
sigmaNMohr = 100 ; 
tauMohr = C0 + muStatic * sigmaNMohr ; 

%   for all directions in 3-space, calculate normal and shear stresses 
%   on all 3d surfaces    
phi_index = 0 ; 
theta_index = 0 ; 
increment = 1 ;
phi_min = 90 ; 
phi_max = 180 ; 
theta_min = 0 ; 
theta_max = 360 ; 
phi_n = (phi_max - phi_min)/increment + 1 ; 
theta_n = (theta_max - theta_min)/increment + 1 ; 
sigmaN = zeros(theta_n, phi_n) ; 
tau = zeros(theta_n, phi_n) ; 
for phi = 90:increment:180

    phi_index = phi_index + 1 ; 
    phi_rad = ( phi - 90 ) * pi / 180 ; 

    theta_index = 0 ; 

    for theta = 0:increment:360

        theta_index = theta_index + 1 ;
        theta_rad = theta * pi / 180 ;
        
        %   convert pole to strike and dip 
        [ strike, dip ] = Pole(theta_rad, phi_rad, 0) ; 
        
        %   calculate normal and shear stress on the plane 
        [ stressFracture, dcStressFracture, ~ ] = ShearOnPlane(stressTensor, trendS1rad, plungeS1rad, trendS3rad, strike, dip) ; 
    
        %   save normal and shear stresses for later calculation 
        sigmaN(theta_index, phi_index) = stressFracture(1,1) ; 
        tau(theta_index, phi_index) = stressFracture(3,1) ; 

    end ; 
    
%   end for each direction 
end ; 

%   calculate tendencies - slip, dilatation and frac. suscep
%   calculate normalised slip tendency (Morris et al., 1996)
TsMax = max(max( tau ./ sigmaN )) ; 
Ts = ( tau ./ sigmaN ) / TsMax ; 

%   calculate dilatation tendency (Ferril et al., 1999)
Td =  ( sorted_sigma(1) - sigmaN ) ./ ( sorted_sigma(1) - sorted_sigma(3) ) ;  

%   calculate shear stress/ dilation displacement ratio 
%   from Delaney et al.(1988)
TD = tau ./ (sigma3 + Pf) ;

%   calculate fracture susceptibility
Sf = sigmaN - ( tau ./ muStatic ) ;    

%   calculate muOA (opening angle), Jolly & Sanderson, 1997
OA = tau ./ (Pf - sigmaN) ;
muOAfracture = atand(OA) ;

%   Calculate and display various stress ratios:
Phi = ( sigma2 - sigma3 ) / ( sigma1 - sigma3 ) ; 
R = (sigma1 - sigma2) / (sigma1 - sigma3) ;
Rprime = ( Pf - sigma3 ) / ( sigma1 - sigma3 ) ; 

disp(' ') ; 
disp(['Stress ratio phi: ', num2str(Phi)]) ; 
disp(['Stress ratio R: ' , num2str(R)]) ;
disp(['Stress ratio R'': ', num2str(Rprime)]) ; 

%   plot azimuthal variation of tendencies, with poles to fractures
%   overlain
deltaP = increment * pi / 180 ;
phiP = pi/2:deltaP:pi ; 
phiP = phiP - pi/2 ; 
thetaP = 0:deltaP:2*pi ; 
[phiP, thetaP] = meshgrid(phiP, thetaP) ;

%   equal area projection 
dp = sqrt(1 - sin(phiP)) ; 
xeqarea = dp .* sin(thetaP) ; 
yeqarea = dp .* cos(thetaP) ;  
rPrim = 1 ; 
xPrim = -rPrim:0.0001:rPrim ;
yPrim = sqrt(rPrim^2 - xPrim.^2) ; 

%  convert Fracture pole plunges and plunge directions to cartesian coords    
dp = sqrt(1 - sin(polesFracturesRad(:,1))) ;
newTrend = rem(polesFracturesRad(:,2), 2*pi) ; 
xFractures = dp .* sin(newTrend) ; 
yFractures = dp .* cos(newTrend) ; 

%   convert principal stress orientations into cartesian coords for equal
%   area
[ pstress, dCp ] = PrincipalStress(stressTensor, trendS1rad, plungeS1rad, trendS3rad) ; 
trendS2rad = pstress(2,2) ; 
plungeS2rad = pstress(2,3) ; 
plungeS3rad = pstress(3,3) ; 
[ xS1, yS1 ] = StCoordLine(trendS1rad, plungeS1rad, 1) ;  
[ xS2, yS2 ] = StCoordLine(trendS2rad, plungeS2rad, 1) ;  
[ xS3, yS3 ] = StCoordLine(trendS3rad, plungeS3rad, 1) ;  

lwPrim = 1 ; 
sizePoleMarker = 15 ; 
sizeStressMarker = 10 ; 
ncontours = 20 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   plot figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   1. stereogram of Slip Tendency, Ts
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

contourf(xeqarea, yeqarea, Ts, ncontours, 'EdgeColor', 'none') ; 
hold on ; 
plot(xPrim, yPrim, '-k', 'LineWidth', lwPrim) ; 
plot(xPrim, -yPrim, '-k', 'LineWidth', lwPrim) ;
plot(xFractures, yFractures, '.r', ...
        'MarkerSize', sizePoleMarker ) ;
plot(xS1, yS1, 's', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
plot(xS2, yS2, 'd', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
plot(xS3, yS3, '^', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
hold off ; 
cmocean(('thermal'), 20) ; 
title({'Slip tendency, T_s';['n=', num2str(nFractures)]}) ; 
view(0, 90) ;  
axis equal off ;
xlim([-1.05 1.05]) ; 
ylim([-1.05 1.05]) ; 
caxis([0 1]) ; 
cb = colorbar ; 
cb.Location = 'SouthOutside' ; 
cb.Label.String = 'Normalised slip tendency' ; 
print -r600 -dtiff 'FracTend_Ts_stereo.tif' ; 

%   2. stereogram of Dilatation Tendency, Td
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

contourf(xeqarea, yeqarea, Td, ncontours, 'EdgeColor', 'none') ; 
hold on ; 
plot(xPrim, yPrim, '-k', 'LineWidth', lwPrim) ; 
plot(xPrim, -yPrim, '-k', 'LineWidth', lwPrim) ;
plot(xFractures, yFractures, '.r', ...
        'MarkerSize', sizePoleMarker ) ;
plot(xS1, yS1, 's', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
plot(xS2, yS2, 'd', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
plot(xS3, yS3, '^', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
hold off ; 
cmocean(('thermal'), 20) ; 
title({'Dilatation tendency, T_d';['n=', num2str(nFractures)]}) ; 
view(0, 90) ;  
axis equal off ;
xlim([-1.05 1.05]) ; 
ylim([-1.05 1.05]) ; 
caxis([0 1]) ; 
cb = colorbar ; 
cb.Location = 'SouthOutside' ; 
cb.Label.String = 'Dilatation tendency' ; 
print -r600 -dtiff 'FracTend_Td_stereo.tif' ; 

%   3. stereogram of fracture susceptibility, Sf
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

contourf(xeqarea, yeqarea, Sf, 24,'EdgeColor', 'none') ; 
hold on ; 
plot(xPrim, yPrim, '-k', 'LineWidth', lwPrim) ; 
plot(xPrim, -yPrim, '-k', 'LineWidth', lwPrim) ;
plot(xFractures, yFractures, '.r', 'MarkerSize', sizePoleMarker ) ; 
plot(xS1, yS1, 's', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
plot(xS2, yS2, 'd', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
plot(xS3, yS3, '^', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
hold off ; 
% Flip the colourbar: most likely to reactivate = warmer colours
colormap(flipud(cmocean(('thermal'), 20))) ;      
title({'Fracture susceptibility, S_f';['n=', num2str(nFractures)]}) ; 
view(0, 90) ;  
axis equal off ;
xlim([-1.05 1.05]) ; 
ylim([-1.05 1.05]) ; 
cb = colorbar ; 
cb.Location = 'SouthOutside' ;  
cb.Label.String = '\DeltaP_{f}, MPa' ; 
print -r600 -dtiff 'FracTend_Sf_stereo.tif' ; 

%   4. stereogram of opening angle, OA
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

cbOA = cmocean('thermal', 18) ; 
cbOA(19, :) = [ 1, 1, 1 ] ;         %   white for OA < 0 
cbOA(20, :) = [ 1, 1, 1 ] ;         %   white for OA < 0 
cbOA = flipud(cbOA) ;               %   reverse it 

contourf(xeqarea, yeqarea, muOAfracture, ncontours, 'EdgeColor', 'none') ; 
hold on ; 
plot(xPrim, yPrim, '-k', 'LineWidth', lwPrim) ; 
plot(xPrim, -yPrim, '-k', 'LineWidth', lwPrim) ;
plot(xFractures, yFractures, '.r', 'MarkerSize', sizePoleMarker) ;
plot(xS1, yS1, 's', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
plot(xS2, yS2, 'd', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
plot(xS3, yS3, '^', ...
        'MarkerSize', sizeStressMarker, 'MarkerEdgeColor','k', 'MarkerFaceColor', 'w') ; 
hold off ; 
colormap(cbOA) ; 
title({['Opening angle \mu_a for P_f=', num2str(Pf), ' MPa']; ['n=', num2str(nFractures)]}) ; 
view(0, 90) ;  
axis equal off ;
xlim([-1.05 1.05]) ; 
ylim([-1.05 1.05]) ; 
caxis([-10 90]) ;
cb = colorbar ; 
cb.Location = 'SouthOutside' ; 
cb.Label.String = 'Opening angle, \circ' ; 
print -r600 -dtiff 'FracTend_OA_stereo.tif' ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   for each plane in the input file
%   new loop to calculate specifc values for the supplied poles 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
sigmaNFracture = zeros(nFractures,1) ; 
tauFracture = zeros(nFractures,1) ; 
for i = 1:nFractures
    
    %   convert fracture pole to strike and dip 
    [ strike, dip ] = Pole(polesFracturesRad(i,2), polesFracturesRad(i,1), 0) ; 
    
    %   calculate normal and shear stress on the plane 
    [ stressFracture, dcStressFracture, R ] = ShearOnPlane(stressTensor, trendS1rad, plungeS1rad, trendS3rad, strike, dip) ; 
    
    %   save normal and shear stresses for later calculation 
    sigmaNFracture(i) = stressFracture(1,1) ; 
    tauFracture(i) = stressFracture(3,1) ; 
    
%   end for each plane 
end ; 

%   calculate normalised slip tendency
TsFractureFile = ( tauFracture ./ sigmaNFracture ) / TsMax ; 
        
%   calculate dilatation tendency
TdFractureFile = ( sorted_sigma(1) - sigmaNFracture ) ./ ( sorted_sigma(1) - sorted_sigma(3) ) ;    

%   calculate fracture susceptibility
SfFractureFile = sigmaNFracture - ( tauFracture ./ muStatic ) ;

%   calculate opening angle
OAFile = tauFracture ./ (Pf - sigmaNFracture) ;
muOAfractureFile = atand(OAFile) ;

%   write out text file of data values for the specific fracture poles
fidData = fopen('PolesWithValues.txt', 'wt') ; 
for i = 1:nFractures 
    fprintf(fidData, '%2.2f, %3.2f, %1.3f, %1.3f, %5.2f, %5.2f, %5.2f, %5.3f\n', ...
                    [ polesFractures(i,1), polesFractures(i,2), ...
                      TsFractureFile(i), TdFractureFile(i), SfFractureFile(i), tauFracture(i), sigmaNFracture(i), muOAfractureFile(i) ] ) ; 
end ; 
fclose(fidData) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  plot Mohr diagrams, with fractures & contoured stability measures  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

sizePoleMarker = 10 ; 

%   1. Ts 
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1) ; 
hold on ; 
contourf(sigmaN, tau, Ts, ncontours, 'EdgeColor', 'none') ; 
plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k') ; 
cmocean('thermal') ; 
for f = 1:nFractures
    plot(sigmaNFracture(f), tauFracture(f), '.r', ...
            'MarkerSize', sizePoleMarker) ; 
end ; 
hold off ; 
xlim([0 sigma1*1.05]) ; 
ylim([0 sigmad*0.75]) ;
xlabel('Effective normal stress, MPa') ; 
ylabel('Shear stress, MPa') ; 
title({'Slip tendency, T_s';['n=', num2str(nFractures)]}) ; 
caxis([0 1]) ; 
cb = colorbar ; 
cb.Location = 'SouthOutside' ; 
cb.Label.String = 'Normalised slip tendency' ; 
print -r600 -dtiff 'FracTend_Ts_mohr.tif' ; 

%   2. Td 
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1) ; 
hold on ; 
contourf(sigmaN, tau, Td, ncontours, 'EdgeColor', 'none') ; 
plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k') ; 
cmocean('thermal') ; 
for f = 1:nFractures
    plot(sigmaNFracture(f), tauFracture(f), '.r', ...
            'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r') ; 
end ; 
hold off ; 
xlim([0 sigma1*1.05]) ; 
ylim([0 sigmad*0.75]) ;
xlabel('Effective normal stress, MPa') ; 
ylabel('Shear stress, MPa') ; 
title({'Dilatation tendency, T_d';['n=', num2str(nFractures)]}) ; 
caxis([0 1]) ; 
cb = colorbar ; 
cb.Location = 'SouthOutside' ; 
cb.Label.String = 'Dilatation tendency' ; 
print -r600 -dtiff 'FracTend_Td_mohr.tif' ; 
 
%   3. Sf
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1) ; 
hold on ; 
contourf(sigmaN, tau, Sf, ncontours, 'EdgeColor', 'none') ; 
plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k') ; 
colormap(flipud(cmocean('thermal'))) ;
for f = 1:nFractures
    plot(sigmaNFracture(f), tauFracture(f), '.r', ...
            'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r') ; 
end ; 
hold off ; 
xlim([0 sigma1*1.05]) ; 
ylim([0 sigmad*0.75]) ;
xlabel('Effective normal stress, MPa') ; 
ylabel('Shear stress, MPa') ; 
title({'Fracture susceptibility, S_f';['n=', num2str(nFractures)]}) ; 
cb = colorbar ; 
cb.Location = 'SouthOutside' ; 
cb.Label.String = '\DeltaP_{f}, MPa' ; 
print -r600 -dtiff 'FracTend_Sf_mohr.tif' ; 

%   opening angle Mohr diagram 
figure ; 
set(gcf, 'PaperPositionMode', 'manual') ; 
set(gcf, 'PaperUnits', 'inches') ; 
set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

%   modified colourbar for opening angle scale 
cbOA = cmocean('thermal', 9) ; 
cbOA(10, :) = [ 1, 1, 1 ] ;         %   white for OA < 0 
cbOA = flipud(cbOA) ;               %   reverse it 

plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1) ; 
hold on ; 
fill([0, Pf, Pf, 0], [0, 0, sigma1, sigma1], 'b', 'FaceAlpha', 0.2) ; 
contourf(sigmaN, tau, muOAfracture, ncontours, 'EdgeColor', 'none') ; 
plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k') ; 
for f = 1:nFractures
    plot(sigmaNFracture(f), tauFracture(f), '.r', ...
            'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r') ; 
end ; 
hold off ; 
colormap(cbOA) ;
xlim([0 sigma1*1.05]) ; 
ylim([0 sigmad*0.75]) ;
xlabel('Effective normal stress, MPa') ; 
ylabel('Shear stress, MPa') ; 
title({['Opening angle \mu_a for P_f=', num2str(Pf), ' MPa']; ['n=', num2str(nFractures)]}) ; 
caxis([-10 90]) ;
cb = colorbar ; 
cb.Location = 'SouthOutside' ; 
cb.Label.String = 'Opening angle, \circ' ; 
print -r600 -dtiff 'FracTend_OA_mohr.tif' ; 
 
disp(' ') ; 
disp(['*** ...finished FracTend version ', num2str(version), ' at ', datestr(now), '.']) ; 
disp(' ') ; 
