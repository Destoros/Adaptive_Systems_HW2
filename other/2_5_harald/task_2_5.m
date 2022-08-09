clear all;
close all;
clc;

load('data.mat');

realization_count=size(X,1);

%task 2.5.a configuration
alg_type_list=[ 0 1 ];%0=Standard LMS, 1=Gradient Descent
mu_list=[0.0001 0.001 0.01 1];
P_list=[size(X,2)];
use_analytic_list=[1];%matters only for gradient descent
N_list=[length(h)];

%task 2.5.b configuration
% alg_type_list=[ 1 ];%0=Standard LMS, 1=Gradient Descent
% mu_list=[ 0.0005 ];
% P_list=[10 100 1000];
% use_analytic_list=[0 1];%matters only for gradient descent
% N_list=[length(h)];

%task 2.5.c configuration
% alg_type_list=[ 0 1 ];%0=Standard LMS, 1=Gradient Descent
% mu_list=[0.001];
% P_list=[size(X,2)];
% use_analytic_list=[1];%matters only for gradient descent
% N_list=[10 100];

for N=N_list
    for use_analytic=use_analytic_list
        for P=P_list
            %only required for task 2.5.a
            cos_x_pos=0:P-1;
            cos_x_pos=pi/2.0*cos_x_pos;
            rxx_analytic_pos=0.5*cos(cos_x_pos);
            rxx_analytic_pos(1,1)=rxx_analytic_pos(1,1)+sigma_u^2;
            rxx_analytic_N=rxx_analytic_pos(1,1:N);
            R_XX_analytic=toeplitz(rxx_analytic_N);%compute analytic R_xx
            p_analytic=R_XX_analytic*h;%compute analytic p

            for alg_type=alg_type_list
                if alg_type==0
                    alg_type_str='Standard LMS';
                else
                    alg_type_str='Gradient Descent';
                end
                for mu=mu_list
                    v={};
                    e={};
                    %if 1==2
                    alpha=0;
                    col_count = 0;
                    for r=1:1:realization_count
                        fprintf("Computing: %s, µ=%f, realization=%d/%d, P=%d, analytic=%d\n",alg_type_str, mu,r,realization_count, P, use_analytic);
                        x_r=X(r,:);
                        d_r=D(r,:);
                        if alg_type==0
                            [y,e_r,c]=lms_algorithm(x_r,d_r,N,mu,alpha,0);
                        else
                            if use_analytic==1
                               R_xx_used=R_XX_analytic;
                               p_used=p_analytic;
                            else
                                %only required for task 2.5.b
                                x_r=X(r,1:P);%take vector of length P (max 10000)
                                d_r=D(r,1:P);%take vector of length P (max 10000)
                                rdx=r_dx(d_r,x_r,P);%compute cross correlation, approximate p
                                rxx=r_dx(x_r,x_r,P);%compute auto correlation, approximate Rxx
                                r_dx_N=rdx(1,P:P+N-1);%cut N values starting at P (P=center of cross corr,"k=0")
                                r_xx_N=rxx(1,P:P+N-1);%cut N values starting at P (P=center of auto corr,"k=0")
                                R_xx_used=toeplitz(r_xx_N);%create auto corr of x
                                p_used=r_dx_N;

                            end
                            [y,e_r,c]=gd_algorithm(x_r,d_r,N,mu,R_xx_used,p_used);
                        end
                        col_count=size(c,2);
                        v_r=[];
                        for col=1:1:col_count
                            c_i=c(:,col);
                            c_i_minus_h=c_i-h; 
                            v_r=[v_r c_i_minus_h];
                        end
                        v{r}=v_r;
                        e{r}=e_r.^2;
                    end%realization_count loop

                    %compute E()
                    E_v=[];
                    E_e=[];
                    for col=1:1:col_count
                        v_avg=zeros(N,1);
                        e_avg=zeros(N,1);
                        for r=1:1:realization_count
                            v_r = v{r};
                            v_r_i=v_r(:,col);
                            v_avg = v_avg+v_r_i;
                            e_r=e{r};
                            e_r_i=e_r(1,col);
                            e_avg = e_avg+e_r_i;
                        end
                        v_avg = 1.0/realization_count*v_avg;
                        e_avg = 1.0/realization_count*e_avg;
                        E_v=[E_v v_avg];
                        E_e=[E_e e_avg];
                    end%col_count loop

                    E_v_0=abs(E_v(:,1));
                    ln_E_v_n_ratio_over_time=[];
                    dbg=[];
                    for col=1:1:col_count
                       E_v_n=abs(E_v(:,col));
                       E_v_n_ratio=E_v_n./E_v_0;
                       %dbg=[dbg E_v_i_ratio];
                       ln_E_v_n_ratio=log(E_v_n_ratio);
                       ln_E_v_n_ratio_over_time=[ln_E_v_n_ratio_over_time ln_E_v_n_ratio];
                    end
                    E_e_0=E_e(:,1);
                    ln_E_e_n_ratio_over_time=[];
                    for col=1:1:col_count
                        E_e_n=E_e(:,col);
                        E_e_n_ratio=E_e_n./E_e_0;
                        ln_E_e_n_ratio=log(E_e_n_ratio);
                        ln_E_e_n_ratio_over_time=[ln_E_e_n_ratio_over_time ln_E_e_n_ratio];
                    end

                    %plot ln(|E(v_k(n))|/|E(v_k(0))|)
                    legend_entries=[];
                    for i=1:1:N
                        legend_entry=sprintf("k=%d",i);
                        legend_entries=[legend_entries legend_entry];
                    end
                    figure;
                    plot(transpose(ln_E_v_n_ratio_over_time));
                    title(sprintf("%s, µ=%f, v_k, P=%d, analytic=%d, N=%d",alg_type_str, mu, P, use_analytic, N))
                    xlabel('n');
                    ylabel('ln(|E(v_k(n))|/|E(v_k(0))|)');
                    legend(legend_entries);

                    %plot ln(E(e^2(n))/E(e^2(0)))
                    figure;
                    plot(transpose(ln_E_e_n_ratio_over_time));
                    title(sprintf("%s, µ=%f, e(n), P=%d, analytic=%d, N=%d",alg_type_str, mu, P, use_analytic, N))
                    xlabel('n');
                    ylabel('ln(E(e^2(n))/E(e^2(0)))');
                    %legend("");

                end%mu loop end
            end %alg_type loop (LMS, GD)
        end%P_list loop
    end%use_analytic list loop
end%N_list loop
%figure(1)
% plot(ln_E_div_E');
% grid
% xlabel('n');
% ylabel('E');
% xv=linspace(0,1,size(ln_E_div_E,2));
% figure(2)
% plot(xv,ln_E_div_E);
% grid
% xlabel('n');
% ylabel('E');
