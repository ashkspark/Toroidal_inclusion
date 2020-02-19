clear all; clc;

global  Ri Ro mu omega Bo


method = 'dqm';

Ri=0.5;        %\delta p/mu=p.
Ro=1.0;
pc = Ro/Ri;
Bo=2.0;
omega=0.01; 
mu=1;
sigma = 1e-6;
delta = 1e-2;
delta1 = delta;
%-----eigenstrains
%{ d omega=1/2*(1-erf((R-Ri)/sigma)); sigma=1e-5;
%}
%--------------------------------------------------------------------------
n = 41; %number of grids in the circumferential direction.
m = 41; %number of grids in the radial direction
%--------------------------------------------------------------------------

switch method
    case 'fdm'
        x = linspace(0,pi,n)';  % = phi
        y = linspace(delta,Ro,m)';  % = R
        
        Dix = FDM(n,x(2)-x(1));
        Diy = FDM(m,y(2)-y(1));
        
        Ix = eye(n);
        Iy = eye(m);
        O    = zeros(m*n);
        D0   = kron(Iy,Ix);
        D1x  = kron(Iy,Dix);
        D1y  = kron(Diy,Ix);
        D2x  = kron(Iy,Dix^2);
        D2y  = kron(Diy^2,Ix);
        D2xy = kron(Diy,Dix);
        
    case 'dqm'
        x = chebspace(0,pi,n)';  % = phi
        y = chebspace(delta,Ro,m)';  % = R
        
        Dix = DQ(x,2);
        Diy = DQ(y,2);
        
        Ix = eye(n);
        Iy = eye(m);
        O    = zeros(m*n);
        D0   = kron(Iy,Ix);
        D1x  = kron(Iy,Dix{1});
        D1y  = kron(Diy{1},Ix);
        D2x  = kron(Iy,Dix{2});
        D2y  = kron(Diy{2},Ix);
        D2xy = kron(Diy{1},Dix{1});  
        
    case 'hdq'
        
        x = chebspace(0,pi,n)';  % = phi
        y = chebspace(delta,Ro,m)';  % = R
        
        Dix = HDQ(x,2);
        Diy = HDQ(y,2);
        
        Ix = eye(n);
        Iy = eye(m);
        O    = zeros(m*n);
        D0   = kron(Iy,Ix);
        D1x  = kron(Iy,Dix{1});
        D1y  = kron(Diy{1},Ix);
        D2x  = kron(Iy,Dix{2});
        D2y  = kron(Diy{2},Ix);
        D2xy = kron(Diy{1},Dix{1});
        
    case 'mixed'
        
        x = chebspace(0,pi,n)';  % = phi
        y = chebspace(delta,Ro,m)';  % = R
        
        Dix = HDQ(x,2);
        Diy = DQ(y,2);
        
        Ix = eye(n);
        Iy = eye(m);
        O    = zeros(m*n);
        D0   = kron(Iy,Ix);
        D1x  = kron(Iy,Dix{1});
        D1y  = kron(Diy{1},Ix);
        D2x  = kron(Iy,Dix{2});
        D2y  = kron(Diy{2},Ix);
        D2xy = kron(Diy{1},Dix{1});   
end

ofs = zeros(m*n,1);

[X,Y] = ndgrid(x,y);

%%

sx1  = Bo+Y(:).*cos(X(:));
sx2  = 3.*Bo+4.*Y(:).*cos(X(:));
%-------------------Equation 3.75
A11 = diag(2./Y(:).^2)-diag(2./Y(:))*D1y-2*D2y-diag(1./Y(:).^2)*D2x...
      -diag(2.*cos(X(:))./sx1)*D1y+diag(1./Y(:).*sin(X(:))./sx1)*D1x...
      +diag(2.*cos(X(:)).^2./sx1.^2);

A12 = -D2xy+diag(2./Y(:))*D1x+diag(Y(:).*sin(X(:))./sx1)*D1y...
      -diag(Y(:).*sin(2.*X(:))./sx1.^2);

A13 = D1y;
%------------------Equation 3.76
A21 = diag(2.*sin(X(:))./Y(:).^2./sx1)-diag(1./Y(:).*sin(2.*X(:))./sx1.^2)...
      -diag(1./Y(:).^3.*sx2./sx1)*D1x-diag(1./Y(:).^2)*D2xy;

A22 = diag(2./Y(:).*sin(X(:))./sx1)*D1x+diag(2.*sin(X(:)).^2./sx1.^2)...
      -diag(1./Y(:).*sx2./sx1)*D1y-D2y-diag(2./Y(:).^2)*D2x;

A23 = diag(1./Y(:).^2)*D1x;
%------------------Equation 3.78
A31 = diag(sx1)+diag(Y(:).*cos(X(:)))+diag(Y(:).*sx1)*D1y;

A32 = -diag(Y(:).^2.*sin(X(:)))+diag(Y(:).*sx1)*D1x;

A33 = O;

A = [A11,A12,A13; A21,A22,A23; A31,A32,A33];  %creating the 3mn*3mn matrix

cn=[+2.*cos(X(:))./sx1.^2;
    -2./Y(:).*sin(X(:))./sx1.^2;
    Y(:)];

F = zeros(3*m*n+1,1);
F(2*m*n+1:3*m*n,1) = +3./2.*Y(:).*sx1.*(1./2.*(1-erf((Y(:)-Ri)./sigma))).*omega;

A = [[A,cn];zeros(1,3*m*n+1)];


% db=0.003;
%     
% F =[-2.*cos(X(:))./sx1.^2.*db;
%     +2./Y(:).*sin(X(:))./sx1.^2.*db;
%     +3./2.*Y(:).*sx1.*(1./2.*(1-erf((Y(:)-Ri)./sigma))).*omega-Y(:).*db]; 

%% BCs
%-----------------------------------
id = reshape(1:m*n,[n,m]);
l1 = id(:,1);    l2 = id(:,2);                                % R = 0;
r1 = id(:,end);  r2 = id(:,end-1);                            % R = R0;
t1 = id(1,:)';   t2 = id(2,:)';                               % phi = 0;     
b1 = id(end,:)'; b2 = id(end-1,:)';                           % phi = pi;
lqid = id((n+1)/2,1);
%-------------boundary condition 3.81b--------------------------------------------

% U
BC = [D0,O,O,ofs];
A(l1+2*m*n,:) = BC(l1,:);
F(l1+2*m*n) = 0;

% BC = [D1y,O,O];
% A(l1+m*n,:) = BC(l1,:);
% F(l1+m*n) = 1./2.*(3.*omega./2-db./Bo);

BC = [D1x,O,O,ofs];      
A(t1,:) = BC(t1,:);
A(b1,:) = BC(b1,:);
F([t1;b1]) = 0;

% BC = [D2xy,O,O];
% A(t2,:) = BC(t1,:);
% A(b2,:) = BC(b1,:);
% F([t2;b2]) = 0;

% BC = [D1x,O,O];
% A(l1,:) = BC(l1,:);
% F(l1) = 0;
%----------------------------------
% W

BC = [O,D0,O,ofs];
A(l1+m*n,:) = BC(l1,:);
F(l1+m*n) = 0;



BC = [O,D0,O,ofs];
A(t1+m*n,:) = BC(t1,:);
A(b1+m*n,:) = BC(b1,:);
F([b1+m*n; t1+m*n]) = 0;



% BC = [O,D2x,O];
% A(t2,:) = BC(t1,:);
% A(b2,:) = BC(b1,:);
% F([t2;b2]) = 0;



% BC = [O,D2x,O];
% A(t2+m*n,:) = BC(t1,:);
% A(b2+m*n,:) = BC(b1,:);
% F([b2+m*n; t2+m*n]) = 0;

%----------------------------
%p

% BC = [O,O,D0];
% A(l1+m*n,:) = BC(l1,:);
% F(l1+m*n) = 2.1393;


% BC = [O,O,D1x];
% A(l1+2*m*n,:) = BC(l1,:);
% F(l1+2*m*n) = 0;


BC = [O,O,D1x,ofs];
A(t1+2*m*n,:) = BC(t1,:);
A(b1+2*m*n,:) = BC(b1,:);
F([b1+2*m*n; t1+2*m*n]) = 0;

%-----------------------------
%mixed

BC = [2*D1y, O, -D0, ofs];
A(r1,:) = BC(r1,:);
F(r1) = 0;


BC = [1./Ro.^2.*D1x, D1y, O, ofs];
A(r1+m*n,:) = BC(r1,:);
F(r1+m*n) = 0;

 
% BC = [diag(2./Bo./delta.^2.*sin(X(:)))-diag(sin(2.*X(:))./delta./Bo.^2)-3./delta.^3.*D1x-1./delta.^2.*D2xy, diag(2./delta./Bo*sin(X(:)))*D1x-3./delta.*D1y-D2y-2./delta.^2.*D2x, 1./delta.^2.*D1x];
%   size(BC)
% A(l1+m*n,:) = BC(l1,:);
% F(l1+m*n) = 2./Bo.^2.*db./delta.*sin(x);


lqu = 2./Bo./delta1.^2.*D0-3./delta1.^3.*D1x-1./delta1.^2.*D2xy;   lqw = 2./Bo./delta1.*D1x-3./delta1.*D1y-D2y-2./delta1.^2.*D2x;   lqp = 1./delta1^2.*D1x;
lqb = -2./Bo.^2./delta1.*ones(m*n,1);
BC = [lqu,lqw,lqp,lqb];
A(3*m*n+1,:) = BC(lqid,:);
F(3*m*n+1)=0;

%% Solution

U = A\F;
db= U(end)

% rank(A)- m*n*3

U = reshape(U(1:3*m*n,1),[],3);

u = U(:,1); %another reshape is required to make the vector a grid.
w = U(:,2);
p = U(:,3);

%% Illustration

figure(1); clf
subplot(1,3,1);
surf(Y.*cos(X), Y.*sin(X),reshape(u,[n,m])); hold on
surf(Y.*cos(X),-Y.*sin(X),reshape(u,[n,m]));

subplot(1,3,2);
surf(Y.*cos(X), Y.*sin(X), reshape(w,[n,m])); hold on
surf(Y.*cos(X),-Y.*sin(X),-reshape(w,[n,m]));

subplot(1,3,3);
surf(Y.*cos(X), Y.*sin(X),reshape(p,[n,m])); hold on
surf(Y.*cos(X),-Y.*sin(X),reshape(p,[n,m]));


