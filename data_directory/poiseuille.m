%% Script to calculate solution to Poiseuille flow at time t
%% THIS FUNCTION IS BASED ON WATERS AND KING 1973
function poiseuille(t,eta,lamda,betaa)
clear u

% flow parameters
ro = 1.0;  % density
%lamda = 1.0;  % relaxation time

theta = betaa*lamda; % retardation time
nu = eta/ro;  
dpdx = 8.0;  % pressure gradient
h = 1.0;  % channel width

% cross channel coord vector, y, and initialise u
y = 0:h/100:h;
u0 = dpdx*h*h/(8.0*eta);

% non-dimensional parameters
y1 = y./h;
t1 = eta*t/(ro*h*h);
S1 = eta*lamda/(ro*h*h);
S2 = eta*theta/(ro*h*h);

u = -4.0.*y1.*(y1.-1);


nn = 30;
for n=1:nn
N=(2*n-1)*pi;

if(lamda>0)
X1 = sin(N*y1)/N**3;
aN=1+S2*N**2;
bN = sqrt((1+S2*N**2)**2 - 4*S1*N**2);
X2 = exp(-(aN*t1)/(2.0*S1));
GN1 = min(cosh(bN*t1/(2.0*S1)),1e300);    %% Avoid GN1 and GN2 going to +Inf
GN2a = (1 + (S2 - 2*S1)*N**2)/bN;
GN2b = min(sinh(bN*t1/(2*S1)),1e300);
GN = GN1 + GN2a*GN2b;

u = u .- 32*X1*X2*GN;

else
X1 = sin(N*y1)/N**3;
X2 = exp(-N*N*t1);
u = u .- 32*X1*X2;
endif

end

u=u.*u0;

plot(y,u,'k-','linewidth',2)


endfunction
