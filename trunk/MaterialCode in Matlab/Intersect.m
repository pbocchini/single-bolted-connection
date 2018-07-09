function [u,f] = Intersect(k1,u1,f1,k2,u2,f2)
%function to find the intersection point between two lines
        u=((f2-k2*u2)-(f1-k1*u1))/(k1-k2);
        f=k1*u+(f1-k1*u1);
        %or f=k2*u+(f2-k2*u2);
end


