clear all;
Z=50;
M=zeros(Z+1,Z+1);
N=6;npg=3;c=0.1;b=1;r=0.3;T=0.2;w=0.8;ro=5;F=1/(1-w);sigma=0.3;
mu=0.02;
s=2;

M(Z+1,Z+1)=1-mu;M(Z+1,Z)=mu;M(1,2)=mu;M(1,1)=1-mu;C=[];

for i=2:Z
    k=i-1;%%%population c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Pd=0;  
  Pc=0;
    for j_c=0:N-1    
        j_d=N-j_c-1;
       
       if j_c>k-1
            E1=0;
       else
            E1=nchoosek(k-1, j_c);
       end
       if (N-j_c-1)>(Z-k)
            E2=0;
       else
            E2=nchoosek(Z-k, N-j_c-1);
       end
       if j_c>k
            E3=0;
       else
            E3=nchoosek(k, j_c);
       end
       if (N-j_c-1)>(Z-k-1)
            E4=0;
       else
            E4=nchoosek(Z-k-1, N-j_c-1);
       end
         if j_c+1<npg
            theta1=0;
        else
            theta1=1;
        end
        if j_c<npg
            theta2=0;
        else
            theta2=1;
        end
      if r>=T && j_c+1<npg
        pi_c=(-c+(1-r)*b)*(ro-1)+(b-c)*(F-ro+1);
      else
        pi_c=(b*theta1+b*(1-r)*(1-theta1)-c)*F;
      end
      if r>=T && j_c<npg
        pi_d=((1-r)*b)*(ro-1)+(b-c)*(F-ro+1)-sigma; 
      else
        pi_d=(b*theta2+b*(1-r)*(1-theta2))*F-sigma;
      end
      
       temp1=E1*E2/nchoosek(Z-1, N-1);
       temp2=E3*E4/nchoosek(Z-1, N-1);
       
        Pc=Pc+temp1*pi_c;
        Pd=Pd+temp2*pi_d;
    end   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_jia=k*(Z-k)/((1+exp(s*(Pd-Pc)))*Z^2);
    T_jian=k*(Z-k)/((1+exp(s*(Pc-Pd)))*Z^2);
    M(i,i-1)=(1-mu)*T_jian+mu*k/Z;
    M(i,i+1)=(1-mu)*T_jia+mu*(Z-k)/Z;
    M(i,i)=1-M(i,i-1)-M(i,i+1);
end
    [x,lumda]=eig(M'); %对其转置即可求得是左特征向量
    r1=abs(sum(lumda));
    q=find(r1==max(r1));
%     max_lumda=lumda(q,q); %最大特征值1
    max_x=x(:,q);%对应的特征向量
    PP=max_x/sum(max_x);
%    x1=0:0.01:1;
plot(PP)


%   [V,D]=eig(M');
%   [x1,x2]=find(D>0.999999);
%   x=V(:,x2);
%   x=x/sum(x);
% 
% plot(x)
% hold on