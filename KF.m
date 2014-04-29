% Filename: KF
% GP fitting with state-based model approach to Kepler Data

function [Xpre,NLL] = KF(Y, phi, sigma, Z, p);
    if nargin < 5;
        p = 6;
    end;

    n = length(Y);

    % Use the Särkkä routines to compute the transition matrices for the filter;
    [F, L, H, Qc, P_inf] = materngp2ss(phi, sigma, p);
    [G, Q] = lti_disc(F, L,Qc);

    % Initialisation
    X = Y(1).*H';
    P = P_inf;
    Xpre = [];
    Ppre = [];
    NLL  =  0.5*log(2*pi);
    for k = 1:n;

        % Predictions w/ t-1 measurements : X(t|t-1), P(t|t-1)
        X = G*X;
        P = G * P * G' + Q;

        % NLL computation
        NLL = NLL + 0.5*log(P(1,1)+Z)+0.5*((Y(k)-X(1)).^2)/(P(1,1)+Z);

        % Kalman gain
        K = P*H'/(P(1,1)+Z);

        % Predictions : X(t|t),P(t|t)
        Xp = X + K * ( Y(k)-X(1));
        Pp = (eye(p+1) - K * H) * P;
        Xpre = [Xpre , Xp(1)];
        Ppre = [Ppre , Pp];

        % Update
        X = Xp;
        P = Pp;
    end
end
