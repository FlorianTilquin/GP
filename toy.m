n = 500;
u = 1:n;

%% Input Parameters;
l = 20;
sigma1 = 2;
sigma2 = 0.5;

%% Output Parameters;
phi = 1/(l*sqrt(2.));
sigma = sigma1;
Z = sigma2;

%% Data points;
% Mean Vector and Covariance matrix;
mean = zeros(1,n);
cov = SE(u,u,sigma1,l,sigma2);

% Drawns from the covariance matrix;
X = mvnrnd(mean, cov,1);
figure();
plot(X,'.k');hold on;

%% Prediction inputs test;
[Xp,~] = KF(X,phi,sigma1,sigma2);

plot(Xp,'-b');
