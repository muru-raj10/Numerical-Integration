function A=GaussQuad2(f,a,b,tol)
x1=-sqrt(1/3);
x2=sqrt(1/3);
%since c1=c2=1, we do not need to define it
c=(a+b)/2;

s1ac=(c-a)/2;       %translation
s2ac=(c+a)/2;
Aac=s1ac*sum(feval(f,s1ac*x1+s2ac)+feval(f,s1ac*x2+s2ac));

s1cb=(b-c)/2;
s2cb=(b+c)/2;
Acb=s1cb*sum(feval(f,s1cb*x1+s2cb)+feval(f,s1cb*x2+s2cb));

s1ab=(b-a)/2;
s2ab=(b+a)/2;
Aab=s1ab*sum(feval(f,s1ab*x1+s2ab)+feval(f,s1ab*x2+s2ab));

Bound = 3*tol;
if abs(Aab-(Aac+Acb))<Bound
    A=Aac+Acb;
else
    A=GaussQuad2(f,a,c,tol/2)+GaussQuad2(f,c,b,tol/2);
end

end