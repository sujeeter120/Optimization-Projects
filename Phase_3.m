n_var = input(" Enter Number of Variable : ");

for i=1:n_var    
x0=zeros(n_var,1);
x0(i)=input('Enter intial guess:');

end
Solution_at_x = Constr_Opt(x0) %calling main function










function x_star=Constr_Opt(x0)
t=1;                              % iteration counter

R=0.1;                             %intial Penality term
e=0.001;                           % Termination tolerance

 x0=phase_2(x0,R)                  %Finding best point at intail guess
 f0=Function(x0,R);
 function_value(t)=f0;

while true
     R=10*R;                       %Updating Penality by facor 10
   x0=phase_2(x0,R);               % Finding optimum point at updated penality     
   f1=Function(x0,R);
   
   E=abs(f1-f0);                   %termination condition for penality method
 
   if (E<e)
       break;
   end
 f0=f1                           % Swapping function value
 t=t+1
 
end
Iteration=[1:t]
  x_star=x0;
 function_value(t)=Function(x0,R)
 plot(Iteration,function_value);
%  xlabel('No of Iteration');
%  ylabel('Function Value');
%  legend('Function Value vs Iteration');
%  

end
 
 
 
function answer =phase_2(x0,R)      % x0 = Initial Approximation 
 P=x0;
 K = 1 ;                                   % Iteration Count
 function_value(K)=Function(x0,R);            %Function Value at intial guess 
  f=1;                                                                    
 n_var= length(x0);
 
 e = 0.001 ;                               %  termination Condition
 
condition= true;
while condition
 
s1=eye(n_var,n_var);                         % N independent search Direction
 x=zeros(n_var,n_var);                        %for storing new point
for j=1:n_var                                 %First loop for geting alpha value
[alfa,F] = Bounding_Search(x0,s1,j,f,R);       % Calling bounding Phase getting alfa value and number of function evalution
  
 for i=1:n_var                               % second loop for getting value of new point by searching all direction for particular  alpha
 
 x(:,j)=x0 + alfa*s1(:,j);                    
 
 
 end
 
 x0=x(:,j);                                   %Updating old point with new Point
end
[m,Fp] = Bounding_Search(x0,s1,1,f,R);            %One more Search in first direction 
                
x1=x0+m*s1(:,1);
X=x1;
Y=x(:,1);
 
d= X-Y;                                        %getting New direction According to Powell_conjugate
p=abs(d'*d);
 q=sqrt(p);
 s2=d/q;                                        % New Search Direction                                      
 LI=(s2'*s1(:,1));                          % Checking of Linear Indepedency
 
 
 
 if(q<e)                                                   %Termination Condition
     break;
 elseif(LI==0)                                   %Condition for checking is our direction is dependent or not. if dependent then terminate loop
     break;
 elseif(K>1000)                                  %No of Iteration
     break;
 else
     for k=n_var:-1:2                           % updating Search direction
     s1(:,k)=s1(:,k-1);
     end
                       % First iteration Completed and we get new point and new search direction.
 end
 s1(:,1)=s2;
 n=Bounding_Search(x0,s1,1,f,R);  
 P=x0+n*s1(:,1);
 answer=P;
 K=K+1;                                 


 
                                                                
end

 answer=P;
%    F=F+Fp;                                            
%  fprintf("The No. of Function evalution:%d\n",F)
%  fprintf("The No. of iterations are : %d \n",K)
%  fprintf("The minimum function value is : %f \n",Function(X))
%  iteration = [1:K]  ;
%  plot(iteration,function_value);
%  xlabel('No of Iteration');
%  ylabel('Function Value');
%  legend('Function Value vs Iteration');
%  
end



% Please Uncomment for Running Code
% ------- Multivariable Function ------------

 function fun_val = Function(x,R)
 %% Question 1
% y = (x(1)-10)^3 +(x(2)-20)^3;
% constr1=(x(1)-5)^2+(x(2)-5)^2-100;
% c1=min(0,constr1);
% constr2=-(x(1)-6)^2-(x(2)-5)^2+82.81;
% c2=min(0,constr2);
% constr3=x(1)-13;
% c3=min(0,constr3);
% constr4=20-x(1);
% c4=min(0,constr4);
% constr5=x(2);
% c5=min(0,constr5);
% constr6=4-x(2);
% c6=min(0,constr6);
% penality= R*(c1^2+c2^2+c3^2+c4^2+c5^2+c6^2);
% fun_val= y +penality;
%% Question 2
y=-((sin(2*pi*x(1)))^3*sin(2*pi*x(2)))/(x(1)^3*(x(1)+x(2)));
constr1=-x(1)^2+x(2)-1;
c1=min(0,constr1);
constr2=x(1)-(x(2)-4)^2-1;
c2=min(0,constr2);
constr3=x(1);
c3=min(0,constr3);
constr4=10-x(1);
c4=min(0,constr4);
constr5=x(2);
c5=min(0,constr5);
constr6=10-x(2);
c6=min(0,constr6);
penality= R*(c1^2+c2^2+c3^2+c4^2+c5^2+c6^2);
fun_val=y + penality;
%% Question 3
%  y=x(1)+x(2)+x(3);
% constr1=1-0.0025*(x(4)+x(6));
% c1=min(0,constr1);
% constr2=1-0.0025*(-x(4)+x(5)+x(7));
% c2=min(0,constr2);
% constr3=1-0.01*(-x(6)+x(8));
% c3=min(0,constr3);
% constr4=-100*x(1)+x(1)*x(6)-833.33252*x(4)+83333.33;
% c4=min(0,constr4);
% constr5=-x(2)*x(4)+x(2)*x(7)+1250*x(4)-1250*x(5);
% c5=min(0,constr5);
% constr6=-x(3)*x(5)+x(3)*x(8)+2500*x(5)-1250000*x(5);
% c6=min(0,constr6);
% constr7=x(1)-100;
% c7=min(0,constr7);
% constr8=10000-x(1);
% c8=min(0,constr8);
% constr9=x(2)-1000;
% c9=min(0,constr9);
% constr10=10000-x(2);
% c10=min(0,constr10);
% constr11=x(3)-1000;
% c11=min(0,constr11);
% constr12=10000-x(3);
% c12=min(0,constr12);
% constr13=x(4)-10;
% c13=min(0,constr13);
% constr14=1000-x(4);
% c14=min(0,constr14);
% constr15=x(5)-10;
% c15=min(0,constr15);
% constr16=1000-x(5);
% c16=min(0,constr16);
% constr17=x(6)-10;
% c17=min(0,constr17);
% constr18=1000-x(6);
% c18=min(0,constr18);
% constr19=x(7)-10;
% c19=min(0,constr19);
% constr20=1000-x(7);
% c20=min(0,constr20);
% constr21=x(8)-10;
% c21=min(0,constr21);
% constr22=1000-x(8);
% c22=min(0,constr22);
% 
% penality=R*(c1^2+c2^2+c3^2+c4^2+c5^2+c6^2+c7^2+c8^2+c9^2+c10^2+c11^2+c12^2+c13^2+c14^2+c15^2+c16^2+c17^2+c18^2+c19^2+c20^2+c21^2+c22^2);
% 
% fun_val=y + penality;





end
 
 
%Single variable Function creation from multivariable
function fun_val=Objective_Fun(y,x0,s1,j,R)     % taking argument y(intial point for bounding),x0( taking variable),s1(direction matrix),j=index)
 
 x0=x0 + y*s1(j,:)'; 
fun_val=Function(x0,R);
end
 
 
 
function [l,F] = Bounding_Search(x0,s1,j,f,R)   % Bounding Phase and Bisection Method for finding alfa
 
 
delta = 0.8;                            % Intial Value taken for finding alfa
 feval = 0;
 y0=0.9;
 
 y2 = y0; y1 = y2 - abs(delta); y3 = y2 + abs(delta);
 f1 = Objective_Fun(y1,x0,s1,j,R);
 f2 = Objective_Fun(y2,x0,s1,j,R);
 f3 = Objective_Fun(y3,x0,s1,j,R);
 i = 1;
 feval= feval+3;
  if(f3<f2 && f2<f1)
     delta = abs(delta);                     % we are checking our increment is positive or negative for bounding
     y1 = y3;
     f1 = f3;
  elseif(f1==f2 && f2==f3)
     delta=0;
 
  else
     delta = -abs(delta);
     
 end
 
while (f1<f2)
    
       f2 = f1;                         % according to sign of delta we swapping f1,f2 and f3
       y3 = y2;
       y2 = y1;
    
       y1 = y2 + (2^i)*delta;              %finding new point according to bounding phase rule
       f1 = Objective_Fun(y1,x0,s1,j,R);
     
       i = i+1;
        feval = feval + 1;
        if(f1>f2)
         break;
        end
end
 
[l,F]=Bisection_search(y1,y3,x0,s1,j,f,R);          % to get better value alfa we aplly bisection method after bounding phase
 
end
 
 
function [alfa,F]=Bisection_search(y1,y3,x0,s1,j,f,R) 
%  fprintf('Bisection Search Method');
 if(y1<y3)
     y1=y1;            % swaping y1 and y3 to get get lower limit and upper limit for bisection method
 else
     t=y1;
     y1=y3;
     y3=t;
 end
 delta_x = .0001; e = 10^-5;           %setting termination condition
 y1  = y1 + delta_x; y2 = y1 -delta_x;    % setting the value to get numerical differentiation for bisection method
 y3 = y3+delta_x; y4 = y3 - delta_x;
 z = (y1 + y3)/2; z1 = z + delta_x; z2 = z - delta_x;
 
 f1 = Objective_Fun(y1,x0,s1,j,R);
 
 f2 = Objective_Fun(y2,x0,s1,j,R);
 
 
 f3 = Objective_Fun(y3,x0,s1,j,R);
 
 
 f4 = Objective_Fun(y4,x0,s1,j,R);
 f5 =  Objective_Fun(z1,x0,s1,j,R);
 f6 =  Objective_Fun(z2,x0,s1,j,R);
  
 
 p = (f1-f2)/(2*delta_x);            %  formula for getting differentiation by central difference formula
 q = (f3-f4)/(2*delta_x);
 r = (f5-f6)/(2*delta_x);
  k =  Objective_Fun(z,x0,s1,j,R);
 feval = 6;
 
 
 
 
 while(abs(r)>=e)                 %checking termination condition
     
 
 if(r>0)
     y1 = y1;
     y3 = z;
    
     
 else
     y1 = z;
     y3 =y3;
     
 end
     
  z = (y1+y3)/2;                               % to get new point
  z1 = z + delta_x; z2 = z - delta_x;
  f5 =  Objective_Fun(z1,x0,s1,j,R);             %updating function for next iteration
  f6 =  Objective_Fun(z2,x0,s1,j,R);
  r = (f5-f6)/(2*delta_x);
  k =  Objective_Fun(z,x0,s1,j,R);
 
 
 
 
  feval= feval+3;
 end
  
alfa=z;                    % here we get alfa that will pass by bounding phase in powell conjugate
 F=feval+f;
end
 
 



 