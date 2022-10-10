% -----Generate Random initial value and check bunch of them at once-----
n_var = input(" Enter Number of Variable : ") ;
a = input("Enter the lower limit of range : "); % Initial point of the range 
b = input("Enter the upper limit of range ");   % Final point of the range



% for i=1:n_var    
% x0=zeros(n_var,1);
% x(i)=input('Enter intial guess:');
% 
% end
% Solution_x = project_phase_2(x0)

 
 % Testing of code by random Number in given range a & b
for i = 1:9
x0 = (b-a).*rand(n_var,1) + a  % generates random points between a and b 
at_x = project_phase_2(x0) 
end
 
 
% x0 = (b-a).*rand(n_val,1) + a  % generates random points between a and b 
% solution = project_phase_2(x0) 
 
 
 
function answer = project_phase_2(x0)    % x0 = Initial Approximation 
 f=0;
 K = 1 ;                                   % Iteration Counter 
 function_value(K)=Function(x0);                                        
 n_var= length(x0);
 
 e = 0.001 ;                               %  termination Condition
 
condition= true;
while condition
 
s1=eye(n_var,n_var);                         % N search Direction
 x=zeros(n_var,n_var);                        %for storing new point
for j=1:n_var                                 %First loop for geting alpha value
[alfa,F] = Bounding_Search(x0,s1,j,f);     
  
 for i=1:n_var                               % second loop for getting value of new point by searching all direction for particular  alpha
 
 x(:,j)=x0 + alfa*s1(:,j);                    
 
 
 end
 
 x0=x(:,j);                                   %Updating old point with new Point
end
[m,Fp] = Bounding_Search(x0,s1,1,f);            %One more Search in first direction 
                 
x1=x0+m*s1(:,1);
X=x1;
Y=x(:,1);
 
d= X-Y;                                        %getting New direction According to Powell_conjugate
p=norm(d'*d);
 q=sqrt(p);
 s2=d/q;                                        %Normalizing Search direction
 
 
 
 
 if(q<e)                                        %Condition for checking is our direction is dependent or not. if dependent then terminate loop
     break;
 elseif(K>1000)                                  %No of Iteration
     break;
 else
     for k=n_var:-1:2                           % updating Search direction
     s1(:,k)=s1(:,k-1);
     end
 s1(:,1)=s2;
                                % First iteration Completed and we get new point and new search direction.
 end
 m = Bounding_Search(x0,s1,1,f);
R=x0+m*s1(:,1);
 K=K+1;                                 
answer=R; 
 
function_value(K)=Function(R);
 
                                                                
end
 
   F=F+Fp;                                            
 fprintf("The No. of Function evalution:%d\n",F)
 fprintf("The No. of iterations are : %d \n",K)
 fprintf("The minimum function value is : %f \n",Function(R))
 iteration = [1:K]  ;
 plot(iteration,function_value);
 xlabel('No of Iteration');
 ylabel('Function Value');
 legend('Function Value vs Iteration');
 
end



% Please Uncomment for Running Code
% ------- Multivariable Function ------------
function fun_val = Function(x)
 
%--------------QUESTION 1 SUM SQUARES FUNCTION------------------
% n_var = length(x) ;
% fun_val = 0 ;
% for i = 1:n_var
%   fun_val = fun_val + i*x(i)^2 ;
% end
% -------------QUESTION 2 ROSENBROCK FUNCTION --------------------- 
% n_var = length(x) ;
% fun_val = 0 ;
% for i = 1:n_var-1
%   fun_val = fun_val + 100*(x(i+1) - x(i)^2)^2 + (x(i) - 1)^2 ;
% end
 
% ------------QUESTION 3 DIXON PRICE FUNCTION ------------------------
% n_var = length(x) ;
% fun_val = (x(1) - 1)^2  ;
% for i = 2:n_var
%   fun_val = fun_val + i*(2*x(i)^2 - x(i-1))^2 ;
% end
% -------------QUESTION 4 TRID FUNCTION -------------------------
n_var = length(x) ;
fun_val = 0  ;
for i = 1:n_var
  fun_val = fun_val + (x(i)-1)^2 ;
end
for i = 2:n_var
  fun_val = fun_val - x(i)*x(i-1) ;
end
% -----------QUESTION 5 ZAKHAROV FUNCTION ------------------------------
% n_var = length(x) ;
% fun_val = (x(1) - 1)^2  ;
% for i = 2:n_var
%   fun_val = fun_val + i*(2*x(i)^2 - x(i-1))^2 ;
% end
% 
% --------------Himmelblau function----------------------------------
%      fun_val = ((x(1))^2 + x(2) - 11)^2 + (x(1) + (x(2))^2 - 7)^2 ;
 
 
end
 
 
%Function creation for single variable
function fun_val=Objective_Fun(y,x0,s1,j)     % taking argument y(intial point for bounding),x0( taking variable),s1(direction matrix),j=index)
 
 x0=x0 + y*s1(j,:)'; 
fun_val=Function(x0);
end
 
 
 
function [l,F] = Bounding_Search(x0,s1,j,f)   % Bounding Phase and Bisection Method for finding alfa
 
 
delta = 0.1;                            % Intial Value taken for finding alfa
 feval = 0;
 y0=0.9;
 
 y2 = y0; y1 = y2 - abs(delta); y3 = y2 + abs(delta);
 f1 = Objective_Fun(y1,x0,s1,j);
 f2 = Objective_Fun(y2,x0,s1,j);
 f3 = Objective_Fun(y3,x0,s1,j);
 i = 1;
 feval= feval+3;
  if(f3<=f2 && f2<=f1)
     delta = abs(delta);                     % we are checking our increment is positive or negative for bounding
     y1 = y3;
     f1 = f3;
 
  else
     delta = -abs(delta);
     
 
     
 
 end
 
while (f1<f2)
    
       f2 = f1;                         % according to sign of delta we swapping f1,f2 and f3
       y3 = y2;
       y2 = y1;
    
       y1 = y2 + (2^i)*delta;              %finding new point according to bounding phase rule
       f1 = Objective_Fun(y1,x0,s1,j);
     
       i = i+1;
        feval = feval + 1;
        if(f1>f2)
         break;
        end
end
 
[l,F]=Bisection_search(y1,y3,x0,s1,j,f);          % to get better value alfa we aplly bisection method after bounding phase
 
end
 
 
function [alfa,F]=Bisection_search(y1,y3,x0,s1,j,f) 
%  fprintf('Bisection Search Method');
 if(y1<y3)
     y1=y1;            % swaping y1 and y3 to get get lower limit and upper limit for bisection method
 else
     t=y1;
     y1=y3;
     y3=t;
 end
 delta_x = .0001; e = 10^-3;           %setting termination condition
 y1  = y1 + delta_x; y2 = y1 -delta_x;    % setting the value to get numerical differentiation for bisection method
 y3 = y3+delta_x; y4 = y3 - delta_x;
 z = (y1 + y3)/2; z1 = z + delta_x; z2 = z - delta_x;
 
 f1 = Objective_Fun(y1,x0,s1,j);
 
 f2 = Objective_Fun(y2,x0,s1,j);
 
 
 f3 = Objective_Fun(y3,x0,s1,j);
 
 
 f4 = Objective_Fun(y4,x0,s1,j);
 f5 =  Objective_Fun(z1,x0,s1,j);
 f6 =  Objective_Fun(z2,x0,s1,j);
  
 
 p = (f1-f2)/(2*delta_x);            %  formula for getting differentiation by central difference formula
 q = (f3-f4)/(2*delta_x);
 r = (f5-f6)/(2*delta_x);
  k =  Objective_Fun(z,x0,s1,j);
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
  f5 =  Objective_Fun(z1,x0,s1,j);             %updating function for next iteration
  f6 =  Objective_Fun(z2,x0,s1,j);
  r = (f5-f6)/(2*delta_x);
  k =  Objective_Fun(z,x0,s1,j);
 
 
 
 
  feval= feval+3;
 end
  
alfa=z;                    % here we get alfa that will pass by bounding phase in powell conjugate
 F=feval;
end
 
 


