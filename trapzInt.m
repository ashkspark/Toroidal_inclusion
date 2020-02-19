function S = trapzInt(R)
R = R(:)';
S = [R 0 0]-[0 0 R];
S(1)   = []; 
S(end) = [];
S(1)   = R(2)-R(1);
S(end) = R(end)-R(end-1);
S = S/2;
end