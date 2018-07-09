function [x,y]=Bisection(u, f, Kb, REP,u_b0, R, Ki, status,f_b0,Kslip)
% Return the intersection of lin (u, f, Ki) and RichardEquation 
tol=0.00001;
if u>0
low=u;   % tension 
high=u+1;
else 
low=u;   % compressionsion 
high=u-1;
end
guess=(low+high)/2;

maxIterations=100000;
numIterations=0;
 while numIterations<maxIterations
     F1=RichardEquation(REP, f_b0, Kslip, R, Ki, guess, u_b0);
     
     if F1<0  %compression
         F2=f+(guess-u)*Kb;
     end
     if F1>0 %tension 
         F2=f+(guess-u)*Kb;
     end
     
     error=abs(F2)-abs(F1);
     if abs(error) <tol
         break;
     end
     if error >0
         high=guess;
     else
         low=guess;
     end
     guess=(low+high)/2;
     numIterations=numIterations+1;
 end
x=guess;
y=f+(x-u)*Kb;

end