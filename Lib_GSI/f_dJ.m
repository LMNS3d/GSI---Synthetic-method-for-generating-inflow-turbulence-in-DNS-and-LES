
%%% Function: Calculate the first derivative of the integrated variance of
%%% velocity divergence to the phases

% Input:
% aRp -- bIm: Rearranged eigenmodes in y
% d_y:    Grid sizes in y
% p_temp: Phases of the eigenmodes

% Output:
% dJ:     First derivative of the integrated variance of velocity divergence
%         to the phases

function dJ = f_dJ(aRp,rRp,bRp,aIp,rIp,bIp,aRm,rRm,bRm,aIm,rIm,bIm,d_y,p_temp)
    
    nky = size(aRp,2);
    p_temp = reshape(p_temp,[nky,2])';

    p_xp = p_temp(1,:);
    p_yp = p_temp(1,:);
    p_zp = p_temp(1,:);
    
    p_xm = p_temp(2,:);
    p_ym = p_temp(2,:);
    p_zm = p_temp(2,:);
    
    %
    c_xp = cos(p_xp);
    s_xp = sin(p_xp);
    c_yp = cos(p_yp);
    s_yp = sin(p_yp);
    c_zp = cos(p_zp);
    s_zp = sin(p_zp);
    
    c_xm = cos(p_xm);
    s_xm = sin(p_xm);
    c_ym = cos(p_ym);
    s_ym = sin(p_ym);
    c_zm = cos(p_zm);
    s_zm = sin(p_zm);
    
    % Elements of the integrated divergence variance J
    Dsp =  sum( (aRp.*c_xp-aIp.*s_xp) + (rRp.*c_zp-rIp.*s_zp) + (bRp.*s_yp+bIp.*c_yp) , 2 );
    Dcp =  sum( (aRp.*s_xp+aIp.*c_xp) + (rRp.*s_zp+rIp.*c_zp) - (bRp.*c_yp-bIp.*s_yp) , 2 );
    Dsm =  sum( (aRm.*c_xm-aIm.*s_xm) - (rRm.*c_zm-rIm.*s_zm) + (bRm.*s_ym+bIm.*c_ym) , 2 );
    Dcm =  sum( (aRm.*s_xm+aIm.*c_xm) - (rRm.*s_zm+rIm.*c_zm) - (bRm.*c_ym-bIm.*s_ym) , 2 );
    
    % Calculate the first derivative with chain rule
    dJdDsp = 2*d_y.*Dsp;
    dJdDcp = 2*d_y.*Dcp;
    dJdDsm = 2*d_y.*Dsm;
    dJdDcm = 2*d_y.*Dcm;
    
    nky = size(aRp,2);
    
    dJ1 = zeros(2,nky);
    
    dJ1(1,:) = sum( dJdDsp.*(-aRp.*s_xp-aIp.*c_xp ) + dJdDcp.*( aRp.*c_xp-aIp.*s_xp ) ,1 ) ...
             + sum( dJdDsp.*( bRp.*c_yp-bIp.*s_yp ) + dJdDcp.*( bRp.*s_yp+bIp.*c_yp ) ,1 ) ...
             + sum( dJdDsp.*(-rRp.*s_zp-rIp.*c_zp ) + dJdDcp.*( rRp.*c_zp-rIp.*s_zp ) ,1 );
    dJ1(2,:) = sum( dJdDsm.*(-aRm.*s_xm-aIm.*c_xm ) + dJdDcm.*( aRm.*c_xm-aIm.*s_xm ) ,1 ) ...
             + sum( dJdDsm.*( bRm.*c_ym-bIm.*s_ym ) + dJdDcm.*( bRm.*s_ym+bIm.*c_ym ) ,1 ) ...
             + sum( dJdDsm.*( rRm.*s_zm+rIm.*c_zm ) + dJdDcm.*(-rRm.*c_zm+rIm.*s_zm ) ,1 );
    
    dJ1 = dJ1';
    
    dJ = dJ1(:);
    
end


