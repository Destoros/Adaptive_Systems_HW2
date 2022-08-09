%Harald Stiegler, 9330054
%Adaptive Systems, UE, Assignment 2, Task 2.4
function [y,e,c] = gd_algorithm(x,d,N,mu,Rxx,p,c0)
if ~exist('c0','var')
    %c0 not passed -> default zero
    c0=zeros(N,1);
else
    c0=c0(:);%enforce to column vector
end
d=transpose(d(:));%enforce to row vector
x=transpose(x(:));
p=p(:);%enforce to column vector
y=[];
c=[];
e=[];
padded=zeros(1,N-1);
x_padded=[padded x];
c_n_minus_1=c0;
for (sliding_window_start=1:1:length(x))
    x_sliding_window = flip(x_padded(1,sliding_window_start:sliding_window_start+N-1));
    y_n=transpose(c_n_minus_1)*transpose(x_sliding_window);
    y=[y y_n];%append to row vector (because x is also one and y has to be same dimension as x)
    d_n=d(1,sliding_window_start);
    e_n=d_n-y_n;
    e=[e e_n];
    dbg1=Rxx*c_n_minus_1;
    dbg2=p-dbg1;
    dbg3=mu*dbg2;
    c_n=c_n_minus_1+dbg3;

    c_n_minus_1 = c_n;%prepare for next loop
    c=[c c_n];
end
end