function runFracTend(sFile, sPath, ... 
                     nS1, nS2, nS3, nPf, ...
                     nS1Trend, nS1Plunge, nS3Trend, ...
                     nMu, nC0, ...
                     fTsStereo, fTdStereo, fSfStereo, fOAStereo, ...
                     fTsMohr, fTdMohr, fSfMohr, fOAMohr)
%   runFracTend.m - script to plot slip and dilation tendency  
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

%   read in poles to specific fractures; tab-delimited text file, 
%   formatted as plunge then trend 
fnFractures = [ sPath, sFile ] ;  
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
sigma1 = nS1 ;      
sigma2 = nS2 ;      
sigma3 = nS3 ;       
sorted_sigma = [ sigma1, sigma2, sigma3 ] ; 
sigmad = sigma1 - sigma3 ; 
disp(['Principal stresses: ', num2str(sorted_sigma), ' in MPa']) ; 

% pore fluid pressure in MPa, for fracture Opening Angle
% calculations 
Pf = nPf ;

%   note the implicit convention: x//s1, y//s2, z//s3
stressTensor = [ sorted_sigma(1), 0, 0 ; ...
                 0, sorted_sigma(2), 0 ; ... 
                 0, 0, sorted_sigma(3) ] ; 

%   read in stress orientation 
%   e.g. for old case of SHmax azimuth of 135 
        % normal fault system: Trend s1=0, Plunge s1=90 - use s3-trend to change orientation of stress field with s1 as the rotation axis
        % thrust fault: Trend of s3 must be >90 from s1; Ps1=0
        % strike-slip: Trend of s3 must be 90 from Trend of s1
trendS1 = nS1Trend ; 
plungeS1 =  nS1Plunge ;    
trendS3 = nS3Trend ;   
trendS1rad = trendS1 * pi / 180 ; 
plungeS1rad = plungeS1 * pi / 180 ; 
trendS3rad = trendS3 * pi / 180 ; 

%   calculate 'missing' stress orientation data
[ cnS1, ceS1, cdS1 ] = SphToCart(trendS1*pi/180, plungeS1*pi/180, 0) ; 
[ cnS3, ceS3, ~ ] = SphToCart(trendS3*pi/180, 0, 0) ; 

%   vS1 dot product vS3 = 0 - they're orthogonal 
cdS3 = - ( cnS1 * cnS3 + ceS1 * ceS3 ) / cdS1 ; 
if isnan(cdS3) || isinf(cdS3)
    cdS3 = 0 ; 
end 
if cdS3 < 0 
    cdS3 = abs(cdS3) ; 
end 

%   vs1 cross product vS3 = vS2 - it's orthogonal to both 
vS2 = cross([cnS1, ceS1, cdS1], [cnS3, ceS3, cdS3]) ;
vS2(isnan(vS2)) = 0 ; 
vS2(isinf(vS2)) = 0 ; 

cnS2 = vS2(1) ; 
ceS2 = vS2(2) ; 
cdS2 = vS2(3) ; 

[ trendS2rad, plungeS2rad ] = CartToSph(cnS2, ceS2, cdS2) ; 
[ ~, plungeS3rad ] = CartToSph(cnS3, ceS3, cdS3) ; 

disp([cnS3, ceS3, cdS3]) ; 

trendS2 = trendS2rad * 180 / pi ; 
plungeS2 = plungeS2rad * 180 / pi ; 
plungeS3 = plungeS3rad * 180 / pi ; 
if abs(trendS2) < 1e-8
    trendS2 = 0.0 ; 
end 
if abs(plungeS2) < 1e-8 
    plungeS2 = 0.0 ; 
end 
if abs(plungeS3) < 1e-8 
    plungeS3 = 0.0 ; 
end 

% disp('Stress orientation:') ; 
% disp(['Sigma1 plunge/trend - ', num2str(plungeS1, '%02.1f'), '/', num2str(trendS1, '%03.1f')]) ;                     
% disp(['Sigma2 plunge/trend - ', num2str(plungeS2, '%02.1f'), '/', num2str(trendS2, '%03.1f')]) ;                     
% disp(['Sigma3 plunge/trend - ', num2str(plungeS3, '%02.1f'), '/', num2str(trendS3, '%03.1f')]) ;                     

%   coefficient of friction & cohesion 
disp(' ') ; 
muStatic = nMu ; 
C0 = nC0 ; 
sigmaNMohr = 100 ; 
tauMohr = C0 + muStatic * sigmaNMohr ; 

%   for all directions in 3-space, calculate normal and shear stresses 
%   on all 3d surfaces    
phi_index = 0 ; 
increment = 0.5 ;
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
        [ stressFracture, ~, ~ ] = ShearOnPlane(stressTensor, trendS1rad, plungeS1rad, trendS3rad, strike, dip) ; 
    
        %   save normal and shear stresses for later calculation 
        sigmaN(theta_index, phi_index) = stressFracture(1,1) ; 
        tau(theta_index, phi_index) = stressFracture(3,1) ; 

    end 
    
%   end for each direction 
end 

%   calculate tendencies - slip, dilation and frac. suscep
%   calculate normalised slip tendency (Morris et al., 1996)
TsMax = max(max( tau ./ sigmaN )) ; 
Ts = ( tau ./ sigmaN ) / TsMax ; 

%   calculate dilation tendency (Ferril et al., 1999)
Td =  ( sorted_sigma(1) - sigmaN ) ./ ( sorted_sigma(1) - sorted_sigma(3) ) ;  

%   calculate shear stress/ dilation displacement ratio 
%   from Delaney et al.(1988)
% TD = tau ./ (sigma3 + Pf) ;

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
[ pstress, ~ ] = PrincipalStress(stressTensor, trendS1rad, plungeS1rad, trendS3rad) ; 
trendS2rad = pstress(2,2) ; 
plungeS2rad = pstress(2,3) ; 
plungeS3rad = pstress(3,3) ; 
[ xS1, yS1 ] = StCoordLine(trendS1rad, plungeS1rad, 1) ;  
[ xS2, yS2 ] = StCoordLine(trendS2rad, plungeS2rad, 1) ;  
[ xS3, yS3 ] = StCoordLine(trendS3rad, plungeS3rad, 1) ;  

disp('Stress orientation:') ; 
disp(['Sigma1 plunge/trend - ', num2str(plungeS1, '%02.1f'), '/', num2str(trendS1, '%03.1f')]) ;                     
disp(['Sigma2 plunge/trend - ', num2str(plungeS2rad*180/pi, '%02.1f'), '/', num2str(trendS2rad*180/pi, '%03.1f')]) ;                     
disp(['Sigma3 plunge/trend - ', num2str(plungeS3rad*180/pi, '%02.1f'), '/', num2str(trendS3, '%03.1f')]) ;                     

lwPrim = 1 ; 
sizePoleMarker = 15 ; 
sizeStressMarker = 10 ; 
ncontours = 20 ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   plot figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fTsStereo 
    %   1. stereogram of Slip Tendency, Ts
    f = figure ; 
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
    title({'Slip tendency, T_s';['n=', num2str(nFractures)]}) ; 
    view(0, 90) ;  
    axis equal off ;
    xlim([-1.05 1.05]) ; 
    ylim([-1.05 1.05]) ; 
    caxis([0 1]) ; 
    cb = colorbar ; 
    cmocean(('thermal'), ncontours) ; 
    cb.Location = 'SouthOutside' ; 
    cb.Label.String = 'Normalised slip tendency' ; 
    guiPrint(f, 'FracTend_Ts_stereo') ; 
end 

if fTdStereo 
    %   2. stereogram of Dilation Tendency, Td
    f = figure ; 
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
    title({'Dilation tendency, T_d';['n=', num2str(nFractures)]}) ; 
    view(0, 90) ;  
    axis equal off ;
    xlim([-1.05 1.05]) ; 
    ylim([-1.05 1.05]) ; 
    caxis([0 1]) ; 
    cb = colorbar ; 
    cmocean(('thermal'), ncontours) ; 
    cb.Location = 'SouthOutside' ; 
    cb.Label.String = 'Dilation tendency' ; 
    guiPrint(f, 'FracTend_Td_stereo') ; 
end 

if fSfStereo 
    %   3. stereogram of fracture susceptibility, Sf
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

    contourf(xeqarea, yeqarea, Sf, ncontours, 'EdgeColor', 'none') ; 
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
    title({'Fracture susceptibility, S_f';['n=', num2str(nFractures)]}) ; 
    view(0, 90) ;  
    axis equal off ;
    xlim([-1.05 1.05]) ; 
    ylim([-1.05 1.05]) ; 
    cb = colorbar ; 
    cmap = cmocean('thermal', ncontours) ;      
    colormap(flipud(cmap)) ; 
    cb.Location = 'SouthOutside' ;  
    cb.Label.String = '\DeltaP_{f}, MPa' ; 
    guiPrint(f, 'FracTend_Sf_stereo') ; 
end 

if fOAStereo 
    %   4. stereogram of opening angle, OA
    f = figure ; 
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
    guiPrint(f, 'FracTend_OA_stereo') ; 
end 

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
    [ stressFracture, ~, ~ ] = ShearOnPlane(stressTensor, trendS1rad, plungeS1rad, trendS3rad, strike, dip) ; 
    
    %   save normal and shear stresses for later calculation 
    sigmaNFracture(i) = stressFracture(1,1) ; 
    tauFracture(i) = stressFracture(3,1) ; 
    
%   end for each plane 
end 

%   calculate normalised slip tendency
TsFractureFile = ( tauFracture ./ sigmaNFracture ) / TsMax ; 
        
%   calculate dilation tendency
TdFractureFile = ( sorted_sigma(1) - sigmaNFracture ) ./ ( sorted_sigma(1) - sorted_sigma(3) ) ;    

%   calculate fracture susceptibility
SfFractureFile = sigmaNFracture - ( tauFracture ./ muStatic ) ;

%   calculate opening angle
OAFile = tauFracture ./ (Pf - sigmaNFracture) ;
muOAfractureFile = atand(OAFile) ;

%   write out text file of data values for the specific fracture poles
fidData = fopen('PolesWithValues.txt', 'wt') ; 

% write model parameters
fprintf(fidData,  'Model Parameters: \n' ) ;
fprintf(fidData,  'sigma1 (MPa), sigma2 (MPa), sigma3 (MPa), Pf (MPa), phi, Rprime\n' ) ; 
fprintf(fidData,  '%3.2f, %3.2f, %3.2f, %3.2f, %1.2f, %1.2f\n', [ sigma1, sigma2, sigma3, Pf, Phi, Rprime ] ) ; 
fprintf(fidData,  ['sigma1 plunge/trend - ', num2str(plungeS1, '%02d'), '/', num2str(trendS1, '%03d'), '\n'] ) ;
fprintf(fidData,  ['sigma2 plunge/trend - ', num2str(plungeS2, '%02d'), '/', num2str(trendS2, '%03d'), '\n'] ) ;
fprintf(fidData,  ['sigma3 plunge/trend - ', num2str(plungeS3, '%02d'), '/', num2str(trendS3, '%03d'), '\n'] ) ;

% write column headers
fprintf(fidData,  'Model Data: \n' ) ;
fprintf(fidData,  'plunge (deg), trend (deg), Ts, Td, Sf (MPa), tau (MPa), sigmaN (MPa), OpeningAngle (deg)\n' ) ; 

for i = 1:nFractures 
    fprintf(fidData, '%2.2f, %3.2f, %1.3f, %1.3f, %5.2f, %5.2f, %5.2f, %5.3f\n', ...
                    [ polesFractures(i,1), polesFractures(i,2), ...
                      TsFractureFile(i), TdFractureFile(i), SfFractureFile(i), tauFracture(i), sigmaNFracture(i), muOAfractureFile(i) ] ) ; 
end
fclose(fidData) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  plot Mohr diagrams, with fractures & contoured stability measures  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

sizePoleMarker = 10 ; 

if fTsMohr 
    %   1. Ts 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

    plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1) ; 
    hold on ; 
    contourf(sigmaN, tau, Ts, ncontours, 'EdgeColor', 'none') ; 
    plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k') ; 
    for ifrac = 1:nFractures
        plot(sigmaNFracture(ifrac), tauFracture(ifrac), '.r', ...
                'MarkerSize', sizePoleMarker) ; 
    end 
    hold off ; 
    xlim([0 sigma1*1.05]) ; 
    ylim([0 sigmad*0.75]) ;
    xlabel('Effective normal stress, MPa') ; 
    ylabel('Shear stress, MPa') ; 
    title({'Slip tendency, T_s';['n=', num2str(nFractures)]}) ; 
    caxis([0 1]) ; 
    cb = colorbar ; 
    cmocean('thermal') ; 
    cb.Location = 'SouthOutside' ; 
    cb.Label.String = 'Normalised slip tendency' ; 
    guiPrint(f, 'FracTend_Ts_mohr') ; 
end 

if fTdMohr 
    %   2. Td 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

    plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1) ; 
    hold on ; 
    contourf(sigmaN, tau, Td, ncontours, 'EdgeColor', 'none') ; 
    plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k') ; 
    for ifrac = 1:nFractures
        plot(sigmaNFracture(ifrac), tauFracture(ifrac), '.r', ...
                'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r') ; 
    end 
    hold off ; 
    xlim([0 sigma1*1.05]) ; 
    ylim([0 sigmad*0.75]) ;
    xlabel('Effective normal stress, MPa') ; 
    ylabel('Shear stress, MPa') ; 
    title({'Dilation tendency, T_d';['n=', num2str(nFractures)]}) ; 
    caxis([0 1]) ; 
    cb = colorbar ; 
    cmocean('thermal') ; 
    cb.Location = 'SouthOutside' ; 
    cb.Label.String = 'Dilation tendency' ; 
    guiPrint(f, 'FracTend_Td_mohr') ; 
end 

if fSfMohr 
    %   3. Sf
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

    plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1) ; 
    hold on ; 
    contourf(sigmaN, tau, Sf, ncontours, 'EdgeColor', 'none') ; 
    plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k') ; 
    for ifrac = 1:nFractures
        plot(sigmaNFracture(ifrac), tauFracture(ifrac), '.r', ...
                'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r') ; 
    end 
    hold off ; 
    xlim([0 sigma1*1.05]) ; 
    ylim([0 sigmad*0.75]) ;
    xlabel('Effective normal stress, MPa') ; 
    ylabel('Shear stress, MPa') ; 
    title({'Fracture susceptibility, S_f';['n=', num2str(nFractures)]}) ; 
    cb = colorbar ; 
    colormap(flipud(cmocean('thermal'))) ;
    cb.Location = 'SouthOutside' ; 
    cb.Label.String = '\DeltaP_{f}, MPa' ; 
    guiPrint(f, 'FracTend_Sf_mohr') ; 
end 

if fOAMohr 
    %   opening angle Mohr diagram 
    f = figure ; 
    set(gcf, 'PaperPositionMode', 'manual') ; 
    set(gcf, 'PaperUnits', 'inches') ; 
    set(gcf, 'PaperPosition', [ 0.25 0.25 5 5]) ; 

    %   modified colourbar for opening angle scale 
    cbOA = cmocean('thermal', 18) ; 
    cbOA(19, :) = [ 1, 1, 1 ] ;         %   white for OA < 0 
    cbOA(20, :) = [ 1, 1, 1 ] ;         %   white for OA < 0 
    cbOA = flipud(cbOA) ;               %   reverse it 

    plot([0, sigmaNMohr], [C0, tauMohr], '-r', 'LineWidth', 1) ; 
    hold on ; 
    cb = colorbar ; 
    colormap(cbOA) ;
    contourf(sigmaN, tau, muOAfracture, 40, 'EdgeColor', 'none') ; 
    plotMohr_OA(sigma1, sigma2, sigma3, 1, '-k') ; 
    fill([0, Pf, Pf, 0], [0, 0, sigma1, sigma1], 'b', 'FaceAlpha', 0.2) ; 
    for ifrac = 1:nFractures
        plot(sigmaNFracture(ifrac), tauFracture(ifrac), '.r', ...
                'MarkerSize', sizePoleMarker, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r') ; 
    end 
    hold off ; 
    xlim([0 sigma1*1.05]) ; 
    ylim([0 sigmad*0.75]) ;
    xlabel('Effective normal stress, MPa') ; 
    ylabel('Shear stress, MPa') ; 
    title({['Opening angle \mu_a for P_f=', num2str(Pf), ' MPa']; ['n=', num2str(nFractures)]}) ; 
    caxis([-10 90]) ;
    cb.Location = 'SouthOutside' ; 
    cb.Label.String = 'Opening angle, \circ' ; 
    guiPrint(f, 'FracTend_OA_mohr') ; 
end 

end 