function D = FDM(n,h)
D = gallery('tridiag',n,-1,0,1);
D(1,1:3)         = [-3 4 -1];
D(end,end-2:end) = [1 -4 3];
D = D/2/h;
end