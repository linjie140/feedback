clc
clear
%close all
c=1;
F=3;N=5;Z=50;
mu=0.01;
beta=2;
T=8;
% p=0.7;
%%%%%%%%%%%%%%%%%%%%%%%%
i_c=1:1:Z-1;
p_num=[];syms p;I=[];
% w_num=[];syms w;I=[];T=1/(1-w);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(i_c)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Pd=0;  
Pc=0;
    for j_c=0:N-1
        for joc=0:i
       if j_c>i_c(i)-1
            E1=0;
       else
            E1=nchoosek(i_c(i)-1, j_c);
       end
       if (N-j_c-1)>(Z-i_c(i))
            E2=0;
       else
            E2=nchoosek(Z-i_c(i), N-j_c-1);
       end
       if j_c>i_c(i)
            E3=0;
       else
            E3=nchoosek(i_c(i), j_c);
       end
       if (N-j_c-1)>(Z-i_c(i)-1)
            E4=0;
       else
            E4=nchoosek(Z-i_c(i)-1, N-j_c-1);
       end
       if joc>j_c
           E5=0;
       else
           E5=nchoosek(j_c, joc);
       end
       temp1=E1*E2*E5*p^joc*(1-p)^(j_c-joc)/nchoosek(Z-1, N-1);
       temp2=E3*E4*E5*p^joc*(1-p)^(j_c-joc)/nchoosek(Z-1, N-1);
        pi_c=F*c*(j_c+1)/N-c+(F*c*(joc+1)/(j_c+1)-c)*(T-1)*(1+p*(N-1))/N;
        pi_d=F*c*j_c/N+(T-1)*F*c*joc*(1+(1-p)*(N-1))/(N*(j_c+1));
        Pc=Pc+temp1*pi_c;
        Pd=Pd+temp2*pi_d;
        end
    end  
   deltaic=i_c(i)/Z*((Z-i_c(i))/(Z-1))*(1/(1+exp(beta*(Pd-Pc)))-1/(1+exp(beta*(Pc-Pd))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    stupid = solve(deltaic==0,p);
   
    y1=double(stupid);%%转成双精度
    y2=y1(abs(imag(y1))<eps(y1));%% 去掉复数根
    yy=y2(y2>=0&y2<=1);%% 0-1之间的实数根
    if isempty(yy)
        continue;
    else
        I(i)=i_c(i);
        aaa=length(yy);%%0-1之间可能出现多个实数根
        for mm=1:aaa
            p_num(i,mm)=yy(mm);%%应该是至多N个根，可以预先看一下每种比例下有几个根，再决定画几组线
        end
    end


end
aaa1=find(p_num(:,1)~=0);%% 第一组根
plot(p_num(aaa1,1),I(aaa1),'*b')
hold on
aaa2=find(p_num(:,2)~=0);%% 第二组根
plot(p_num(aaa2,2),I(aaa2),'*b')
hold on