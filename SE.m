function K = SE(X,Y,a,b,c);
if nargin < 5;
   c = 0;
end;

XX=sum(X.*X,1);
YY=sum(Y.*Y,1);
XY=X'*Y;
d = abs(repmat(XX',[1 size(YY,2)]) + repmat(YY,[size(XX,2) 1]) - 2*XY);

K = a*exp(-d./(2.*b^2)) + c*eye(length(X));
