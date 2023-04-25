%c-dc
clear all
Z=50;
N=6;M=3;c=0.1;b=1;r=0.3;T=0.2;w=0.8;ro=3;F=1/(1-w);sigma=0.3;
mu=0.02;
beta=2;
i_c=0:1:Z;

deltaic=[];
for i=1:size(i_c,2)
  Pd=0;  
  Pc=0;
    for j_c=0:N-1    
       
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
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if j_c+1<M
            theta1=0;
        else
            theta1=1;
        end
        if j_c<M
            theta2=0;
        else
            theta2=1;
        end
        
      if r>=T && j_c+1<M
        pi_c=((1-r)*b-c)*(ro-1)+(b-c)*(F-ro+1);
      else
        pi_c=(b*theta1+b*(1-r)*(1-theta1)-c)*F;
      end
      
      if r>=T && j_c<M
        pi_d=((1-r)*b)*(ro-1)+(b-c)*(F-ro+1)-sigma; 
      else
        pi_d=(b*theta2+b*(1-r)*(1-theta2))*F-sigma;
      end

       temp1=E1*E2/nchoosek(Z-1, N-1);
       temp2=E3*E4/nchoosek(Z-1, N-1);
       
        Pc=Pc+temp1*pi_c;
        Pd=Pd+temp2*pi_d;
    end   
    deltaic(i)=(1-mu)*(i_c(i)/Z)*((Z-i_c(i))/(Z))*(1/(1+exp(beta*(Pd-Pc)))-1/(1+exp(beta*(Pc-Pd))))+mu*((Z-2*i_c(i))/Z);
%    deltaic(i)=i_c(i)/Z*((Z-i_c(i))/(Z-1))*(1/(1+exp(beta*(Pd-Pc)))-1/(1+exp(beta*(Pc-Pd))));
%    deltaic(i)=(i_e(i)/Z)*((Z-i_e(i))/Z)*(1/(1+exp(beta*(Pd-Pe)))-1/(1+exp(beta*(Pe-Pd))));
%    deltaic(i)=(i_c(i)/Z)*((Z-i_c(i))/Z)*tanh(beta/2*(Pc-Pd));
end

    plot(i_c/Z,deltaic,'*')
    hold on