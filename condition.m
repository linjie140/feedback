function xdot=condition(t,x)
xdot = zeros(2,1); % the number of variable
%N=6;b=1;c=0.1;r=0.3;M=3;T=0.5;w=0.8;ro=1;F=1/(1-w);sigma=0.1;%%FIGURE1%\ro
N=6;b=1;c=0.1;r=0.5;M=3;T=0.2;w=0.8;ro=5;F=1/(1-w);sigma=0.1;%figure2%M
Pi_c=0;
Pi_d=0;
Pi_i=0;
   for Nc=0:N-1
     for Nd=0:N-1-Nc
           Ni=N-Nc-Nd-1;
           temp=factorial(N-1)/(factorial(Nc)*factorial(Nd)*factorial(Ni))*x(1)^Ni*x(2)^Nc*(1-x(1)-x(2))^Nd;  
           if Nc+Ni+1<M
               theta1=0;
           else
               theta1=1;
           end
           if Nc+1<M
               theta2=0;
           else
               theta2=1;
           end
           if Nc+Ni<M
               theta3=0;
           else
               theta3=1;
           end
           if Nc<M
               theta4=0;
           else
               theta4=1;
           end
               if Nc+1<M && r>=T
                      pic = (-c+(1-r)*b)*(ro-1)+(b*theta1-c+(1-r)*b*(1-theta1))*(F-ro+1);
               else
                      pic=(-c+b*theta2+(1-r)*b*(1-theta2))*F;
               end
                  if Nc<M && r>=T
                      pii=(1-r)*b*(ro-1)+(b*theta1-c+(1-r)*b*(1-theta1))*(F-ro+1)-sigma;
                      pid=(1-r)*b*(ro-1)+(b*theta3+(1-r)*b*(1-theta3))*(F-ro+1);
                  else
                      pii=(b*theta4+(1-r)*b*(1-theta4))*F-sigma;
                      pid=(b*theta4+(1-r)*b*(1-theta4))*F;
                  end
          

          Pi_c=Pi_c+temp*pic;
          Pi_d=Pi_d+temp*pid;
          Pi_i=Pi_i+temp*pii;
            
     end
   end

P_=x(1)*Pi_i+x(2)*Pi_c+(1-x(1)-x(2))*Pi_d;
%%%%%%%%%%%%%%%%%%%%%%%%无突变
xdot(1)=x(1)*(Pi_i-P_);
xdot(2)=x(2)*(Pi_c-P_);

%%%%%%%%%%%%%%%%%%%%%%%%有突变
% xdot(1)=x(1)*Pi_c*(1-2*mu)+x(2)*Pi_d*mu+(1-x(1)-x(2))*Pi_i*mu-x(1)*P_;
% xdot(2)=x(2)*Pi_d*(1-2*mu)+x(1)*Pi_c*mu+(1-x(1)-x(2))*Pi_i*mu-x(2)*P_;
