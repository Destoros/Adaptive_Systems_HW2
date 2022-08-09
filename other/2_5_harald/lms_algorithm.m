%Harald Stiegler, 9330054
%Adaptive Systems, UE, Assignment 2, Task 2.4
function [y,e,c] = lms_algorithm(x,d,N,mu,alpha,OPTS,c0)
if ~exist('c0','var')
    %c0 not passed -> default zero
    c0=zeros(N,1);
else
    c0=c0(:);%enforce to column vector
end
d=transpose(d(:));%enforce to row vector, hopatatschig?
x=transpose(x(:));
y=[];
c=[];
e=[];
padded=zeros(1,N-1);
x_padded=[padded x];%pad N-1 zeros at begin for startup: at startup sliding window (length=N) will overlap only by last element with first element of x
c_n_minus_1=c0;
for (sliding_window_start=1:1:length(x))
    x_sliding_window = flip(x_padded(1,sliding_window_start:sliding_window_start+N-1));
    y_n=transpose(c_n_minus_1)*transpose(x_sliding_window);%y=c^T*x
    y=[y y_n];%append to row vector (because x is also one and y has to be same dimension as x)
    d_n=d(1,sliding_window_start);
    e_n=d_n-y_n;
    e=[e e_n];
    update_step = mu*e_n;
    if (OPTS==1)
        update_step=mu*e_n/(norm(x_sliding_window)^2+alpha);
    end
    c_n = c_n_minus_1+transpose(update_step * x_sliding_window);
    c_n_minus_1 = c_n;%prepare for next loop
    c=[c c_n];
    
end

end