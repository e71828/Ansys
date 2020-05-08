clc
clear

GaussP=[-0.774597 0 0.774597];                              %???
h1=1;h2=2;                                           %????
l = 0:h1:h1;
m = 0:h2:h2;
l(2,:)=circshift(l,[0 -1]);
x1 = GaussP(1)*h1/2+sum(l)/2;
x2 = GaussP(2)*h1/2+sum(l)/2;
x3 = GaussP(3)*h1/2+sum(l)/2;

m(2,:)=circshift(m,[0 -1]);
y1 = GaussP(1)*h2/2+sum(m)/2;
y2 = GaussP(2)*h2/2+sum(m)/2;
y3 = GaussP(3)*h2/2+sum(m)/2;


% f =@(x,y) 3*y^2+2*x^3+x*y+x^2*y;
f =@(x,y) 5*y^4;
s = 0;
for i=1:length(x1)-1
    for j=1:length(y1)-1
    s = s+ 5/9*5/9*h1*h2*f(x1(i),y1(j))/4+ 8/9*8/9*h1*h2*f(x2(i),y2(j))/4+8/9*5/9*h1*h2*f(x1(i),y2(j))/4+ 8/9*5/9*h1*h2*f(x2(i),y1(j))/4+...
        5/9*5/9*h1*h2*f(x1(i),y3(j))/4+ 5/9*5/9*h1*h2*f(x3(i),y1(j))/4+8/9*5/9*h1*h2*f(x2(i),y3(j))/4+ 8/9*5/9*h1*h2*f(x3(i),y2(j))/4+...
        5/9*5/9*h1*h2*f(x3(i),y3(j))/4;
    end
end
s