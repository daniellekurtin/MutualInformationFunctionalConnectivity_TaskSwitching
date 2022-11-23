function phi = calcphi(ths) 
% ths must be N x time
N = size(ths,1);
phi = (sum(exp(1i*ths),1))/N;
end