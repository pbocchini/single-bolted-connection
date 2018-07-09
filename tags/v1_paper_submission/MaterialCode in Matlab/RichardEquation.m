function F=RichardEquation(REP, f_b0, Kslip, R, Ki, u, u_b0)
%Compute the stress based on given Richard Equation 

def_norm=abs((u- u_b0)*Ki/(R-abs(f_b0)));  %positive
F=(def_norm*REP(1)/(1+(def_norm*REP(1)/REP(3))^REP(4))^(1/REP(4))+def_norm*REP(2))*(R-abs(f_b0))+abs(f_b0);
if f_b0>0
F=F+Kslip*(u- u_b0);   %positive for tension 
else 
F=-F+Kslip*(u- u_b0);  %negative for compression 
end 
end


