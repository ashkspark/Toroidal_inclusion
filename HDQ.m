function D = HDQ(x,r)

x=x(:)';
nx = numel(x);    a = x(1);     b = x(end);
c=1;
k = c*pi/(b-a);   x = k*(x - a);

ext = ones(1,nx);
DX  = (x'*ext-ext'*x)/2;
I = eye(nx);

q = prod( (sin(DX)+I), 2);

D1 = 0.5*(q*ext)./(sin(DX).*(ext'*q')+I);
D1 = D1 - diag(sum(D1,2));
D1 = D1 - diag(sum(D1,2));

D2 = D1.*( 2*diag(D1)*ext-cot(DX+I) );
D2 = D2 - diag(sum(D2,2));
D2 = D2 - diag(sum(D2,2));

D3 = D1.*(0.5+3*diag(D2)*ext)-1.5*D2.*cot(DX+I);
D3 = D3 - diag(sum(D3,2));
D3 = D3 - diag(sum(D3,2));

D{1}=k*D1;
D{2}=k^2*D2;
D{3}=k^3*D3;

for i = 2: ceil(r/2)
    D{2*i}     = D{i}*D{i};  
    D{2*i+1}   = D{i}*D{i+1};
end

D{end+1} = speye(nx);

end
