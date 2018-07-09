function [f_trial,status_trial,dir_trial,u_transition_trial,f_transition_trial,u_bc_trial,f_bc_trial,u_bt_trial,f_bt_trial,u_b0t_trial,f_b0t_trial,u_b0c_trial,f_b0c_trial]...
    =CalcForceStatusDir(ubkp, fbkp, kbkp,u_trial,status,dir,u_current,f_current,u_st, f_st, u_sc, f_sc, Kslip, Kti, Kci,u_transition,f_transition,...
    u_bc,f_bc,u_bt,f_bt,REPt,REPc,Rt,Rc,u_b0t,f_b0t,u_b0c,f_b0c,Keb,Pslip)
 u_transition_trial=u_transition;
 f_transition_trial=f_transition;
 u_bc_trial=u_bc;
 f_bc_trial=f_bc;
 u_bt_trial=u_bt;
 f_bt_trial=f_bt;
 u_b0t_trial=u_b0t;
 f_b0t_trial=f_b0t;
 u_b0c_trial=u_b0c;
 f_b0c_trial=f_b0c;
 if u_trial<ubkp(1)
     f_trial=RichardEquation(REPc, f_b0c, Kslip, Rc, Kci, u_trial, u_b0c);
 elseif u_trial>=ubkp(1) && u_trial<ubkp(2)
     f_trial=fbkp(2)+kbkp(1)*(u_trial-ubkp(2));
 elseif u_trial>=ubkp(2) && u_trial<ubkp(3)
     f_trial=fbkp(3)+kbkp(2)*(u_trial-ubkp(3));
 elseif u_trial>=ubkp(3) && u_trial<ubkp(4)
     f_trial=fbkp(4)+kbkp(3)*(u_trial-ubkp(4));
 elseif u_trial>=ubkp(4) && u_trial<ubkp(5)
     f_trial=fbkp(5)+kbkp(4)*(u_trial-ubkp(5));    
 elseif u_trial>=ubkp(5) && u_trial<ubkp(6)
     f_trial=fbkp(6)+kbkp(5)*(u_trial-ubkp(6));  
 else
     f_trial=RichardEquation(REPt, f_b0t, Kslip, Rt, Kti, u_trial, u_b0t);
 end

 switch status
     case 0 % current point not in bearing stage 
         if u_trial>ubkp(5) 
             status_trial=1;
         elseif  u_trial<ubkp(2)
             status_trial=-1;     
         else 
             status_trial=status;  %status not change
         end
         
     case 1 % current point in the bearing nonstick tension
         if u_trial>=ubkp(5)
             status_trial=1;
         elseif u_trial<ubkp(5) && u_trial>=ubkp(4)
             status_trial=4;
             if dir==1
                 u_transition_trial=ubkp(5);
                 f_transition_trial=fbkp(5);
             elseif dir==-1
                 u_transition_trial=ubkp(4);
                 f_transition_trial=fbkp(4);
             end
         elseif u_trial<ubkp(4) && u_trial>=ubkp(3)
             status_trial=1;
         elseif u_trial<ubkp(3) && u_trial>=ubkp(2)
             status_trial=0;
             % update the new tension bearing point
             u_bt_trial=Intersect(Keb,ubkp(3),fbkp(3),Kslip,u_st, f_st);
             f_bt_trial=f_st+Kslip*(u_bt_trial-u_st);
             %comporession bearing point upadate !!!!!!!!!!!!!!!!!!!!!!!
             u_bc_trial=(u_bt_trial-u_bt)*0.2+u_bc;   %20% is not going to elongation
             f_bc_trial=f_sc+Kslip*(u_bc_trial-u_sc);
             u_b0c_trial=(u_bt_trial-u_bt)*0.2+u_b0c;
             f_b0c_trial=f_sc+Kslip*(u_b0c_trial-u_sc);
         elseif u_trial<ubkp(2)
             status_trial=-1;
         end
         
     case 4 % current point in the bearing stick tension
         if u_trial>=ubkp(5)
             status_trial=1;
         elseif u_trial<ubkp(5) && u_trial>=ubkp(4)
             status_trial=4;
         elseif u_trial<ubkp(4) && u_trial>=ubkp(3)
             status_trial=1;
         elseif u_trial<ubkp(3) && u_trial>=ubkp(2)
             status_trial=0;
         elseif u_trial<ubkp(2)
             status_trial=-1;
         end
         
     case -1 % current point in the bearing nonstick compression
         if u_trial>=ubkp(5)
             status_trial=1;
         elseif u_trial<ubkp(5) && u_trial>=ubkp(4)
             status_trial=0;
             % update the new compression bearing point 
             u_bc_trial=Intersect(Keb,ubkp(3),fbkp(3),Kslip,u_sc, f_sc);
             f_bc_trial=f_sc+Kslip*(u_bc_trial-u_sc);
             %tension bearing point upadate !!!!!!!!!!!!!!!!!!!!!!!
             u_bt_trial=(u_bc_trial-u_bc)*0.4+u_bt;   %40% is not going to elongation
             f_bt_trial=f_st+Kslip*(u_bt_trial-u_st);
             u_b0t_trial=(u_bc_trial-u_bc)*0.4+u_b0t;
             f_b0t_trial=f_st+Kslip*(u_b0t_trial-u_st);
         elseif u_trial<ubkp(4) && u_trial>=ubkp(3)
             status_trial=-1;
         elseif u_trial<ubkp(3) && u_trial>=ubkp(2)
             status_trial=-4;
             if dir==1
                 u_transition_trial=ubkp(2);
                 f_transition_trial=fbkp(2);
             elseif dir==-1
                 u_transition_trial=ubkp(3);
                 f_transition_trial=fbkp(3);
             end
         elseif u_trial<ubkp(2)
             status_trial=-1;
         end
         
     case -4 % current point in the bearing stick tension
         if u_trial>=ubkp(5)
             status_trial=1;
         elseif u_trial<ubkp(5) && u_trial>=ubkp(4)
             status_trial=0;
         elseif u_trial<ubkp(4) && u_trial>=ubkp(3)
             status_trial=-1;
         elseif u_trial<ubkp(3) && u_trial>=ubkp(2)
             status_trial=-4;
         elseif u_trial<ubkp(2)
             status_trial=-1;
         end   
         
         
            
 end
 if (u_trial>=u_bt_trial && f_trial>0) || (u_trial>=u_bt_trial-2*Pslip/Keb && f_trial<0)   %in tension bearing 
    if f_trial-f_current>=0
        dir_trial=1;
    else
        dir_trial=-1;
    end
 elseif (u_trial<=u_bc_trial && f_trial<0) || (u_trial<=u_bc_trial+2*Pslip/Keb && f_trial>0) %in compression bearing 
    if f_trial-f_current>=0
        dir_trial=-1;
    else
        dir_trial=1;
    end
 else 
     dir_trial=0;
 end 
 

     
     
     
 end


