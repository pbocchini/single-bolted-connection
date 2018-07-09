
function u=history2(n,Upeak)
%history II
%Create Deformation according Upeak 


N=(size(Upeak,2)-1)*n;
u_history2=zeros(1,N+1);
for i=1:size(Upeak,2)-1
    distance=Upeak(i+1)-Upeak(i);
    step=distance/n;
    
    Nstart=(i-1)*n+1;
    Nend=Nstart+n-1;
    
    u_history2(Nstart:Nend)=(Upeak(i):step:Upeak(i+1)-step);
end
u_history2(N+1)=Upeak(end);
u=u_history2;

end



