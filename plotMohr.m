function [ xm, ym ] = plotMohr(s1, s2, s3, lw, sSpec) 

thetaMohr = 0:1:360 ; 
sin2thetaMohr = sind(2 * thetaMohr) ; 
cos2thetaMohr = cosd(2 * thetaMohr) ; 

tau13Mohr = ( ( s1 - s3 ) / 2 ) * sin2thetaMohr ; 
sigma13Mohr = ( s1 + s3 ) / 2 ... 
                       + ( ( s1 - s3 ) / 2 ) * cos2thetaMohr ;  

tau12Mohr = ( ( s1 - s2 ) / 2 ) * sin2thetaMohr ; 
sigma12Mohr = ( s1 + s2 ) / 2 ... 
                       + ( ( s1 - s2 ) / 2 ) * cos2thetaMohr ;  

tau23Mohr = ( ( s2 - s3 ) / 2 ) * sin2thetaMohr ; 
sigma23Mohr = ( s2 + s3 ) / 2 ... 
                       + ( ( s2 - s3 ) / 2 ) * cos2thetaMohr ;  

% lw = 1 ; 
plot(sigma13Mohr, tau13Mohr, sSpec, 'LineWidth', lw) ;
plot(sigma12Mohr, tau12Mohr, sSpec, 'LineWidth', lw) ;
plot(sigma23Mohr, tau23Mohr, sSpec, 'LineWidth', lw) ;

axis equal on ; 
box on ; 
grid on ; 

xm = max(sigma13Mohr) ; 
ym = max(tau13Mohr) ; 
