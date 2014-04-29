% MATERNGP2SS State space representation of GP with Matern kernel
%
% Syntax:
%   [F,L,H,q] = materngp2ss(phi,sigma,p)
% 
% In:
%     phi - Precision parameter (range)
%   sigma - Amount of variation (partial sill)
%       p - Integer part of the smoothness parameter nu = p + 1/2
% Out:
%   F - Feedback matrix
%   L - Noise gain
%   H - Measurement matrix
%   q - Spectral density of noise
%   P_inf - Stationary covariance
%
% Description:
%
%   Form a state space representation to Gaussian process with
%   Matern covariance function
%
%     C(t,t') = sigma^2 exp(-sqrt(2*(p+1)/2) tau/l) Gamma(p+1)/Gamma(2*p+1)
%             * sum_{i=0}^p [(p+i)!/(i!(p-i)!)(sqrt(8(p+1/2))tau/l)^(p-i)]
% 
%   where tau = |t-t'| and l is the length-scale (phi = sqrt(2)/2*l).
%

% Copyright (C) 2010 Jouni Hartikainen, Simo S�rkk�
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function [F,L,H,q,P_inf] = materngp2ss(phi,sigma,p)
    la = 2*sqrt((p+0.5))*phi;
    c = sigma^2.*2*pi^(0.5)*exp(gammaln(p+1)-gammaln(p+0.5))*la.^(2*p+1);
    q = c;

    ppoly = pascal(p+2);
    ppoly = diag(flipud(ppoly))';
    ppoly = ppoly(2:end);
    
    lav = zeros(1,p+1);
    for i = 1:length(lav)
        lav(i) = la^i;
    end

    F = diag(ones(p,1),1);
    F(end,:) = fliplr(-lav.*ppoly); 
    
    L = zeros(p+1,1);
    L(end) = 1;
    H = zeros(1,p+1);
    H(1) = 1;
    P_inf = are(F',zeros(size(F)),L*q*L');
end

