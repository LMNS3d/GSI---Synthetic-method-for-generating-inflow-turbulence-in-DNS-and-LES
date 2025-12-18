
%%% Function: Optimize the phases of eigenmodes with the BFGS algorithm 
%%% such that the integrated variance of velocity divergence is minimized
%%% Application scope: Flow field that is homogeneous in z while
%%% inhomogeneous in y

% Input:
% eig_Up -- eig_Wm: Eigenmodes in y
% ft:       Frequency corresponding to temporal direction
% kz:       Wavenumber corresponding to spatial z direction
% U:        Mean streamwise velocity
% Dy:       First differencing matrix in y
% d_y:      Grid sizes in y
% p_init:   Initialized phases (randomly distributed)

% Output:
% p_opt:    Optimized phases of the eigenmodes
% J:        Integrated variance of velocity divergence corresponding to p_opt (after minimization)
% J_init:   Integrated variance of velocity divergence corresponding to p_init (before minimization)
% rec:      Record of variations of J during the optimization process

function [p_opt,J,J_init,rec] = f_GSI_1D(eig_Up,eig_Vp,eig_Wp,eig_Um,eig_Vm,eig_Wm,ft,kz,U,Dy,d_y,p_init)

    [~,nky] = size(eig_Up);

    p_temp1 = p_init;
    rec = [];
    
    % Divide into two groups regarding sin(ft) and cos(ft)
    kx = -ft./U;
    alppRp = real(kx .* eig_Up);
    gampRp = real(kz .* eig_Wp);
    betpRp = real(Dy  * eig_Vp);
    
    alppIp = imag(kx .* eig_Up);
    gampIp = imag(kz .* eig_Wp);
    betpIp = imag(Dy  * eig_Vp);
    
    alppRm = real(kx .* eig_Um);
    gampRm = real(kz .* eig_Wm);
    betpRm = real(Dy  * eig_Vm);
    
    alppIm = imag(kx .* eig_Um);
    gampIm = imag(kz .* eig_Wm);
    betpIm = imag(Dy  * eig_Vm);
    
    J_init = f_J(alppRp,gampRp,betpRp,alppIp,gampIp,betpIp,alppRm,gampRm,betpRm,alppIm,gampIm,betpIm,d_y,p_temp1);
    
    B = eye(2*nky);
    
    p_temp2 = p_temp1;
    
    count = 1;
    
    J = J_init; vstp = 1;
    
    while J >= 1e-8

        J  = f_J(alppRp,gampRp,betpRp,alppIp,gampIp,betpIp,alppRm,gampRm,betpRm,alppIm,gampIm,betpIm,d_y,p_temp1);
        dJ = f_dJ(alppRp,gampRp,betpRp,alppIp,gampIp,betpIp,alppRm,gampRm,betpRm,alppIm,gampIm,betpIm,d_y,p_temp1);
        
        rec = [rec;J];
        
        %%
        stp = -B\dJ;
        
        %% BFGS
        count2 = 1;
        
        if count <= 50
            
            % Fix the step size in first 50 steps
            
            siz = 0.5;
            p_temp2  = p_temp1 + siz*stp;
            J2  = f_J(alppRp,gampRp,betpRp,alppIp,gampIp,betpIp,alppRm,gampRm,betpRm,alppIm,gampIm,betpIm,d_y,p_temp2);
            dJ2 = f_dJ(alppRp,gampRp,betpRp,alppIp,gampIp,betpIp,alppRm,gampRm,betpRm,alppIm,gampIm,betpIm,d_y,p_temp2);

            p_temp1 = p_temp2;
            
        else
            
            % Determine the step size with Wolfe condition from the 51st step

            c1 = 1e-4;
            c2 = 0.9;
            siz = 100;
            
            while 1

                siz = 0.5*siz;

                p_temp2  = p_temp1 + siz*stp;

                J2  = f_J(alppRp,gampRp,betpRp,alppIp,gampIp,betpIp,alppRm,gampRm,betpRm,alppIm,gampIm,betpIm,d_y,p_temp2);
                dJ2 = f_dJ(alppRp,gampRp,betpRp,alppIp,gampIp,betpIp,alppRm,gampRm,betpRm,alppIm,gampIm,betpIm,d_y,p_temp2);

                phi0  = J;
                dphi0 = stp' * dJ;

                phi   = J2;
                dphi  = stp' * dJ2;

                cd1 = (phi <= phi0 + c1*siz*dphi0);            % Armijo condition
                cd2 = (abs(dphi) <= c2 * abs(dphi0));          % Strong Wolfe condition
                
                count2 = count2 + 1;

                if (cd1 == 1 && cd2 == 1) || (J2 < J) || (count2 >= 1e2) || (siz <= 1e-12)
                    break
                end

            end

            p_temp1 = p_temp2;

            %%% BFGS
            s = siz*stp;
            y = dJ2 - dJ;
            B = B + (y*y')/(y'*s) - (B*(s*s')*B')/(s'*B*s);
            if isnan(trace(B))
                B = eye(2*nky);
            end
            B = B + 1e-2*eye(2*nky);
            
        end
        
        J = J2;
        
        count = count+1;
        if count >= 2e3 || count2 >= 1e2 || siz <= 1e-12
            break
        end

        if (count >= 500 && vstp <= 1e-8)
            break
        end
        
        vstp  = max(abs(stp));

    end
    
    f_J(alppRp,gampRp,betpRp,alppIp,gampIp,betpIp,alppRm,gampRm,betpRm,alppIm,gampIm,betpIm,d_y,p_temp1);
    
    p_opt = reshape(p_temp1,[nky,2])';
    
end

