
%%% Function: Calculate the differencing matrix for first derivative
%%% The code implements strategies from
%%% suggested by W. Don and S. Solomonoff in SIAM J. Sci. Comp.
%%% Vol. 6, pp. 1253--1268 (1994).

% Input:
% yn:   Coordinates
% odr1: Accuracy order
% BC_1: Boundary condition at the starting point of yn
% BC_2: Boundary condition at the end point of yn

% Output:
% Dy:   Differencing matrix for first derivative

function Dy = f_Diff(yn,odr1,BC_1,BC_2)
    
    ny = length(yn);
    
    ell = 1;
    
    Dy = zeros(ny);
    
    for iiy = 1:ny
        
        odr = odr1+2;
       
        ord_m1 = floor(odr/2);
        ord_p1 = ceil(odr/2);

        Dp0 = zeros(1,ny+odr-2 );
        
        yn2 = yn;
        yy2 = [2*yn2(1)-flipud(yn2(2:ord_m1));yn2;2*yn2(end)-flipud(yn2(end-ord_p1+1:end-1))];
        
        % Construct the differencing coefficients based on Taylor series
        H = (yy2(iiy:iiy-1+odr-1)-yy2(iiy+ord_m1-1))/max(abs(yy2(iiy:iiy-1+odr-1)-yy2(iiy+ord_m1-1)));
        mH = zeros(odr-1);
        for iii = 1:odr-1
            mH(iii,:) = H.^(iii-1)/factorial(max(iii-1,1));
        end
        
        tgt = zeros(odr-1,1);
        tgt(ell+1) = 1;
        coeff = (mH\tgt)'/max(abs(yy2(iiy:iiy-1+odr-1)-yy2(iiy+ord_m1-1))).^(ell);
        Dp0(1,iiy:iiy+odr-2 ) = coeff;

        Dy(iiy,:) = Dp0(1,1+(ord_m1-1):end-(ord_p1-1));
        
        % Implement boundary conditions based on BC_1 and BC_2
        
        % Boundary condition at the starting point
        if strcmp(BC_1, 'periodic') == 1      

            Dy(iiy,end-(ord_p1)+2:end) = Dy(iiy,end-(ord_p1)+2:end) + (Dp0(1,1:(ord_p1-1)));

        elseif strcmp(BC_1, 'symmetry') == 1 || strcmp(BC_1, 'no-slip wall') == 1

            Dy(iiy,2:ord_m1) = Dy(iiy,2:ord_m1)-fliplr(Dp0(1,1:(ord_m1-1)));
            
        elseif strcmp(BC_1, 'free') == 1

            Dy(iiy,1)   = Dy(iiy,1) + 2*sum(Dp0(1,1:(ord_m1-1)),2);
            Dy(iiy,2:(ord_m1)) = Dy(iiy,2:(ord_m1))-fliplr(Dp0(1,1:(ord_m1-1)));

        end
        
        % Boundary condition at the end point
        if strcmp(BC_2, 'periodic') == 1

            Dy(iiy,1:(ord_m1-1)) = Dy(iiy,1:(ord_m1)-1) + (Dp0(1,end-(ord_m1-2):end));

        elseif strcmp(BC_2, 'symmetry') == 1 || strcmp(BC_2, 'no-slip wall') == 1

            Dy(iiy,end-(ord_p1)+1:end-1) = Dy(iiy,end-(ord_p1)+1:end-1)-fliplr(Dp0(1,end-(ord_p1-2):end));
            
        elseif strcmp(BC_2, 'free') == 1

            Dy(iiy,end) = Dy(iiy,end) + 2*sum(Dp0(1,end-(ord_p1-2):end),2);
            Dy(iiy,end-(ord_p1)+1:end-1) = Dy(iiy,end-(ord_p1)+1:end-1)-fliplr(Dp0(1,end-(ord_p1-2):end));

        end

    end

end
