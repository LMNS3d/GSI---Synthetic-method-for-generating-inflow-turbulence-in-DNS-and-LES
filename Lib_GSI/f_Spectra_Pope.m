
%%% Function: Model spectrum presented in Turbulent Flows (Stephen B. Pope,
%%% 2000) at Section 6.5.3.

% Input:
% K: Turbulence kinetic energy ((<u'u'>+<v'v'>+<w'w'>)/2)
% L: Length scale
% nu: Kinetic viscosity
% kn: Radius wavenumber

% Output:
% E: Homogeneous isotropic spectrum

function E = f_Spectra_Pope(K,L,nu,kn)

    varepsilon = K^(3/2)./L;
    eta = (nu^3/varepsilon)^(1/4);

    C = 1.5;
    
    p0 = 2;
    beta = 5.2;
    
    c_L = 6.78;
    c_eta = 0.401685;
    
    fl = ((kn.*L)./sqrt((kn.*L).^2 + c_L)).^(5/3 + p0);

    fn = exp(-beta*(((kn.*eta).^4 + c_eta^4).^(1/4) - c_eta));

    E = C*varepsilon.^(2/3).*kn.^(-5/3).*fl.*fn;

end


