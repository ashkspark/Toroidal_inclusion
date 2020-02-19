function x = chebspace(a,b,n)
%(c) Jun-2017 by A Golgoon
if nargin==1, n=a;a=0;b=1; end
x = a+(b-a)*(1-cos(((1:n)-1)*pi/(n-1)))/2;
 end
