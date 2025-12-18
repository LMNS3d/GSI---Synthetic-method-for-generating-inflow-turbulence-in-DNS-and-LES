
%%% Function: Calculate the integrated variance of velocity divergence as 
%%% cost function of the optimization problem

% Input:
% aRp -- bIm: Rearranged eigenmodes in y
% d_y:    Grid sizes in y
% p_temp: Phases of the eigenmodes

% Output:
% J:      Integrated variance of velocity divergence

function J = f_J(aRp,rRp,bRp,aIp,rIp,bIp,aRm,rRm,bRm,aIm,rIm,bIm,d_y,p_temp)

    nky = size(aRp,2);
    p_temp = reshape(p_temp,[nky,2])';

    p_xp = p_temp(1,:);    % 1,nky
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
    Dsp =  sum( (aRp.*c_xp-aIp.*s_xp) + (rRp.*c_zp-rIp.*s_zp) + (bRp.*s_yp+bIp.*c_yp) , 2 );  % ny,nky -> ny,1
    Dcp =  sum( (aRp.*s_xp+aIp.*c_xp) + (rRp.*s_zp+rIp.*c_zp) - (bRp.*c_yp-bIp.*s_yp) , 2 );
    Dsm =  sum( (aRm.*c_xm-aIm.*s_xm) - (rRm.*c_zm-rIm.*s_zm) + (bRm.*s_ym+bIm.*c_ym) , 2 );
    Dcm =  sum( (aRm.*s_xm+aIm.*c_xm) - (rRm.*s_zm+rIm.*c_zm) - (bRm.*c_ym-bIm.*s_ym) , 2 );
    
    J = sum( (Dsp.^2 + Dcp.^2 + Dsm.^2 + Dcm.^2).*d_y , 1 );
    
end
