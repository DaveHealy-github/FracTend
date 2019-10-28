function plotTensor2(x, y, z, t, scale, type, incr)

eigt = eig(t) ; 
a = max(eigt) ; 
b = eigt(2) ; 
c = min(eigt) ;

%   plot tensor t as variety of different glyphs 
it = 0 ; 
for theta = 0:incr:360 
    it = it + 1 ; 
    ip = 0 ; 
    for phi = 0:incr:180 
        ip = ip + 1 ; 
        l(it, ip) = cosd(theta)*sind(phi) ; 
        m(it, ip) = sind(theta)*sind(phi) ; 
        n(it, ip) = cosd(phi) ; 
        %   quadric calculation 
        t1 = t(1,1)*l(it, ip) + t(1,2)*m(it, ip) + t(1,3)*n(it, ip) ; 
        t2 = t(2,1)*l(it, ip) + t(2,2)*m(it, ip) + t(2,3)*n(it, ip) ; 
        t3 = t(3,1)*l(it, ip) + t(3,2)*m(it, ip) + t(3,3)*n(it, ip) ; 
        tn(it, ip) = t1*l(it, ip) + t2*m(it, ip) + t3*n(it, ip) ; 
        trt = ( 1 / abs(t(1,1)) + 1 / abs(t(2,2)) + 1 / abs(t(3,3))) ; 
        tq(it, ip) = sqrt(trt/abs(tn(it, ip))) ;
        if tq(it, ip) > 3.2 * trt 
            tq(it, ip) = 3.2 * trt ; 
        end 
        %   Reynolds 
        tr(it, ip) = tn(it, ip) ; 
        %   HWY
        th(it, ip) = sqrt(t1^2 + t2^2 + t3^2 - tn(it, ip)^2) ;  
        %   PNS
        ct2 = cosd(theta)^2 ; 
        st2 = sind(theta)^2 ; 
        cp2 = cosd(phi)^2 ; 
        sp2 = sind(phi)^2 ; 
        tp(it, ip) = sqrt( 1 / ( ( sp2*ct2 ) / (a^2) + ( sp2*st2 ) / (b^2) + cp2/(c^2) ) ) ; 
    end 
end 

switch type
    case 'Q'
        %   quadric 
        tq = scale .* tq ; 
        surf(x+(l.*tq), y+(m.*tq), z+(n.*tq), tq, 'EdgeColor', 'none') ; 
    case 'R' 
        %   Reynolds, for 'normal' components  
        tr = scale .* tr ; 
        surf(x+(l.*tr), y+(m.*tr), z+(n.*tr), tr, 'EdgeColor', 'none') ; 
    case 'H' 
        %   HWY (after Hashash, Wotring & Yao), for 'shear' components 
        th = scale .* th ; 
        surf(x+(l.*th), y+(m.*th), z+(n.*th), th, 'EdgeColor', 'none') ; 
    case 'P'
        %  PNS (Principal Normal Surface), 
        tp = scale .* tp ; 
        surf(x+(l.*tp), y+(m.*tp), z+(n.*tp), tp, 'EdgeColor', 'none') ; 
    otherwise
        error('') ; 
end 

end 