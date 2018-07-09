function [ubkp, fbkp, kbkp]=CalcBreakpoints(dir, status, Pslip,  u_bc, f_bc, u_bt, f_bt, u_b0c,f_b0c,u_b0t,f_b0t,u_s0c, f_s0c, u_s0t, f_s0t, u_current,f_current, Keb, Kti, Kci, REPc, REPt, Rc, Rt,Kslip, u_transition, f_transition)
tol=10^-7;
ubkp=zeros(6,1);
fbkp=zeros(6,1);
kbkp=zeros(4,1);
%%%%%%%%%%%%%%%%Case 1 u_current before bearing %%%%%%%%%%%%%%%%%%%%%%%%%
 if f_current<=Line(u_current, Keb, u_bc, f_bc) && f_current>=Line(u_current, Keb, u_bt, f_bt)
    if u_bc==u_b0c %there is no  enlongation in the compression yet
          ubkp(1)=u_bc;
          fbkp(1)=f_bc;
    else                            
          [ubkp(1),fbkp(1)]=Bisection(u_bc, f_bc, Kci, REPc,u_b0c, Rc, Kci, status,f_b0c,Kslip); %u_bc intersection with RE at Compression
    end   
    ubkp(2)=u_bc;
    fbkp(2)=f_bc;
    [ubkp(3),fbkp(3)] = Intersect(Kslip,u_s0c,f_s0c,Keb,u_current,f_current);  %intersect between compression slipping and u_current
    [ubkp(4),fbkp(4)] = Intersect(Kslip,u_s0t,f_s0t,Keb,u_current,f_current);  %intersect between compression slipping and u_current
    ubkp(5)= u_bt;
    fbkp(5)=f_bt;
    if u_bt==u_b0t %there is no  enlongation in the tension yet
          ubkp(6)=u_bt;
          fbkp(6)=f_bt;
    else 
          [ubkp(6),fbkp(6)]=Bisection(u_bt, f_bt, Kti, REPt,u_b0t, Rt, Kti, status,f_b0t,Kslip); %u_bt intersection with RE at tension
    end  
    kbkp=[Kci; Kslip; Keb; Kslip; Kti];



%%%%%%%%%%%%%%%%Case 2 u_current after tension bearing %%%%%%%%%%%%%%%%%%%%
 elseif status>0  
    
              if u_bc==u_b0c %there is no  enlongation in the compression yet
                  ubkp(1)=u_bc;
                  fbkp(1)=f_bc;
              else 
                  [ubkp(1),fbkp(1)]=Bisection(u_bc, f_bc, Kci, REPc,u_b0c, Rc, Kci, status,f_b0c, Kslip); %u_bc intersection with RE at Compression
              end   
              ubkp(2)=u_bc;
              fbkp(2)=f_bc;
              if status==4   %stick
                if f_transition>=f_current
                    ubkp(5)=u_transition;            % equals to transition point
                    fbkp(5)=f_transition;
                    ubkp(4)=u_transition-2*Pslip/Keb;             % transition point 2Pslip downward with Keb 
                    fbkp(4)=f_transition-2*Pslip;
                else
                    ubkp(4)=u_transition;            % equals to transition point
                    fbkp(4)=f_transition;
                    ubkp(5)=u_transition+2*Pslip/Keb;             % transition point 2Pslip upward with Keb
                    fbkp(5)=f_transition+2*Pslip; 
                end
   
              elseif status~=4 %non stick
                 if dir ==1   %loading
                  ubkp(4)=u_current-2*Pslip/Keb;             % transition point 2Pslip downward with Keb
                  fbkp(4)=f_current-2*Pslip;   
                  ubkp(5)=u_current;            % equals to current point
                  fbkp(5)=f_current; 
                 elseif dir==-1
                  ubkp(4)=u_current;            % equals to transition point
                  fbkp(4)=f_current; 
                  ubkp(5)=u_current+2*Pslip/Keb;             % transition point 2Pslip upward with Keb
                  fbkp(5)=f_current+2*Pslip;   
                 end 
              end
              if abs(fbkp(5)-RichardEquation(REPt, f_b0t, Kslip, Rt, Kti, ubkp(5), u_b0t))<tol  %break point 5 equals and larger than break point 6
                  ubkp(6)=ubkp(5);
                  fbkp(6)=fbkp(5);
              else          
              [ubkp(6),fbkp(6)]=Bisection(ubkp(5), fbkp(5), Kti, REPt,u_b0t, Rt, Kti, status,f_b0t,Kslip); %u_bt intersection with RE at tension
              end
              [ubkp(3),fbkp(3)] = Intersect(Kslip,u_s0c,f_s0c,Kti,ubkp(4),fbkp(4));  %intersect between compression slipping and point4
              
              kbkp=[Kci; Kslip; Kti; Keb; Kti];
                  
            
              


%%%%%%%%%%%%%%%%Case 3 u_current after compression bearing %%%%%%%%%%%%%%%%    
  elseif status<0  
      if u_bt==u_b0t %there is no  enlongation in the tension yet
                  ubkp(6)=u_bt;
                  fbkp(6)=f_bt;
      else 
                  [ubkp(6),fbkp(6)]=Bisection(u_bt, f_bt, Kti, REPt,u_b0t, Rt, Kti, status,f_b0t, Kslip); %u_bc intersection with RE at tension
      end   
      ubkp(5)=u_bt;
      fbkp(5)=f_bt;
        
      if status==-4  
         if f_transition>=f_current
                    ubkp(3)=u_transition;            % equals to transition point
                    fbkp(3)=f_transition;
                    ubkp(2)=u_transition-2*Pslip/Keb;             % transition point 2Pslip downward with Keb 
                    fbkp(2)=f_transition-2*Pslip;
         else
                    ubkp(2)=u_transition;            % equals to transition point
                    fbkp(2)=f_transition;
                    ubkp(3)=u_transition+2*Pslip/Keb;             % transition point 2Pslip upward with Keb
                    fbkp(3)=f_transition+2*Pslip; 
         end
      
      elseif status~=-4  
          if dir==1  %loading
                ubkp(2)=u_current;            % equals to current point
                fbkp(2)=f_current; 
                ubkp(3)=u_current+2*Pslip/Keb;             % current point 2Pslip downward with Keb
                fbkp(3)=f_current+2*Pslip;
          elseif dir==-1 %unloading
                ubkp(3)=u_current;            % equals to current point
                fbkp(3)=f_current; 
                ubkp(2)=u_current-2*Pslip/Keb;             % transition point 2Pslip downward with Keb
                fbkp(2)=f_current-2*Pslip;
          end
      end
      if abs(fbkp(2)-RichardEquation(REPc, f_b0c, Kslip, Rc, Kci, ubkp(2), u_b0c))<tol %break point 2 equals and smaller than break point 1
                  ubkp(1)=ubkp(2);
                  fbkp(1)=fbkp(2);
      else
                  
                  [ubkp(1),fbkp(1)]=Bisection(ubkp(2), fbkp(2), Kci, REPc,u_b0c, Rc, Kci, status,f_b0c,Kslip); %u_bc intersection with RE at compression
      end
      [ubkp(4),fbkp(4)] = Intersect(Kslip,u_s0t,f_s0t,Kci,ubkp(3),fbkp(3));  %intersect between tension slipping and point 
      
      kbkp=[Kci; Keb; Kci; Kslip; Kti];
  

 end

 
 
 
 
 

end










