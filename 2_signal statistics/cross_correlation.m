function [rdx, mxx] = cross_correlation(x,y,P,N)


if N > P
    error('samples to average P must be greater or equal to filter coefficients N')
end

x_pad = [x(:); zeros(N-1,1)];
y = y(:); %make sure its a col vector

P_window = 1:P;

for k = 0:N-1
    rdx(k+1) = x_pad(P_window + k).' * y(P_window) / P;
end

mxx = 0:N-1;
rdx = rdx(:); %make sure its a col vector


end

% [rdx, mxx] = cross_correlation([1 2 3 4 5],[10 20 30 40 50],3, 2)