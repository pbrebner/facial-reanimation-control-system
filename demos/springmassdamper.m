function dydt=springmassdamper(t,y,m,b,k,f)
global i

dydt=zeros(2,1);
dydt(1)=y(2);
dydt(2)=-(b/m)*y(2)-(k/m)*y(1) + f(i)/m;

end