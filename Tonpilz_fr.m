clear
clc
format long
rhoe=7500;ve=2922;le=0.06;de=0.015;
rhob=8600;vb=sqrt(10.4e10/rhob);
lb=0.06;db=0.06;
rhof=2790;vf=sqrt(7.15e10/rhof);
lf=0.06;df=0.03;dft=0.06;

mileold = 0;
% syms fr
for fr = 1000:0.1:20000
omega=2*pi*fr;
a=3*sqrt(3)/8;%pi/4   1/2 
% a=1/2;
% a=pi/4;
Se=pi/4*de^2;Sb=pi/4*  db^2;Sf=a*df^2;
F=df/(dft-df);
phi0=asin(1/sqrt((omega/vf)^2*lf^2*(F+1)^2+1));
t1=rhoe*ve*Se/(rhob*vb*Sb)*cot(omega/vb*lb);
t2=rhoe*ve*Se/(rhof*vf*Sf)/(tan(omega/vf*lf+phi0)-1/(F*omega/vf*lf));
% t2=rhoe*ve*Se/(rhof*vf*Sf)*cot(omega/vf*lf);
% (t1+t2)/(1-t1*t2)
if abs(tan(omega/ve*le)-(t1+t2)/(1-t1*t2))<0.0005
    fr
    tan(omega/ve*le)-(t1+t2)/(1-t1*t2)
    le1 = atan(t1)*ve/omega
    break;
else
    mile = atan((t1+t2)/(1-t1*t2))*ve/omega ;
    if mile > mileold
        mileold = mile;
        frold = fr;
        le1 = atan(t1)*ve/omega;
    end
end
end
% mileold
% frold
% le1
% vpasolve(tan(omega/ve*le)-(t1+t2)/(1-t1*t2),fr,[7000;20000])%to solve the omega