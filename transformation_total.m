%%%%%%%total transformation%%%%%%%
clear
clc
format short

%%%%%%%coordiate transformation%%%%%%%
% a = eye(3);  %%polariaztion: z

a = [0 1 0;  %%polariaztion: y
     0 0 1
     1 0 0]; 

 
% a = [0 0 1;  %%polariaztion: x
%      1 0 0
%      0 1 0]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%%%%%%%stress transformation%%%%%%%
T(1:6,1:6) =0;
for i=1:3
    for j=1:3
        if i==j 
            alpha = j;
        else
            alpha = 9-i-j;
        end
        for p=1:3
            for q=1:3
                if p==q
                    beta=p;
                else
                    beta =9-p-q;
                end
                if alpha<=3 && beta<=3
                    T(alpha,beta) = a(i,p)*a(i,p);
                end
                if alpha>3 && beta<=3
                    T(alpha,beta) = a(i,p)*a(j,p);
                end
                if alpha<=3 && beta>3
                    T(alpha,beta) = a(i,q)*a(i,p)+a(i,p)*a(j,q);
                end
                if alpha>3 && beta>3
                    T(alpha,beta) = a(i,p)*a(j,q)+a(i,q)*a(j,p);
                end
            end
        end
    end
end
T   ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%%%%%%%strain transformation%%%%%%%
R = diag([1 1 1 2 2 2]);
Tbar = R*T/R ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


%%%%%check%%%%%
check_value = inv(T)/Tbar';
%%%%%%%%%%%%%%%

%%%%%%%total transformation%%%%%%%
syms C c11 c12 c13 c14 c15 c16 c21 c22 c23 c24 c25 c26 c31 c32 c33 c34 c35 c36
syms c41 c42 c43 c44 c45 c46 c51 c52 c53 c54 c55 c56 c61 c62 c63 c64 c65 c66
C = [c11 c12 c13 c14 0 0% c15 c16;
     c21 c22 c23 c24 0 0% c25 c26;
     c31 c32 c33 c34 0 0% c35 c36;
     c41 c42 c43 c44 0 0% c45 c46;
     0 0 0 0 0 0;
     0 0 0 0 0 0;
%      c51 c52 c53 c54 c55 c56;
%      c61 c62 c63 c64 c65 c66;
     ];
C_new = T*C*T';
RT = eye(6);
RT([4,5,6],:) = RT([6,4,5],:); %%ANSYS

C_ANSYS = RT*C_new*RT'
%%%%%the stiffness matrix%%%%%
%%c = RT*T*c*T'*RT'  


%%%%%the compliance matrix%%%%%
%%s = inv(C_ANSYS) = RT*Tbar*s*Tbar'*RT'


%%%%%the piezoelectric strain matrix%%%%%
%%d = a*d*Tbar'*RT'
syms d d15 d31 d33
d = [0   0   0   0   d15 0;
     0   0   0   d15 0   0;
     d31 d31 d33 0   0   0];
d = a*d*Tbar'*RT'; %(3*6)
d_ANSYS = d' %(6*3)

%%%%%the piez oelectric stress matrix%%%%%
%%e = a*e*T'*RT'
syms e e11 e12 e13 e14 e15 e16 e21 e22 e23 e24 e25 e26 e31 e32 e33 e34 e35 e36
e = [e11 e12 e13 e14 0 0;
     e21 e22 e23 e24 0 0;
     0 0 0 0 0 0];
e = a*e*T'*RT'; %(3*6)
e_ANSYS = e' %(6*3)
     

%%%%%the dielectric matrix%%%%%
%%epsilon = a*epsilon*a'
syms epsilon ep11 ep22 ep33
epsilon = diag([ep11 ep22 ep33]);
epsilon = a*epsilon*a'


%%%%%%at constant electric field
%d(3*6) = e*s 
%e = d*c

% cde=[13.2  7.3  7.1  0   0   0;
%       7.3 11.5  7.3  0   0   0;
%       7.1  7.3 13.2  0   0   0;
%       0    0    0    2.6 0   0;
%       0    0    0    0   2.6 0;
%       0    0    0    0   0   3];
% inv(cde)

cde=[13.9  7.43  7.78  0   0   0;
      7.43 11.5  7.43  0   0   0;
      7.78  7.43 13.9  0   0   0;
      0    0    0    2.56 0   0;
      0    0    0    0   2.56 0;
      0    0    0    0   0   3.06]*1e10;
% inv(cde)
% sqrt(1e10/0.1273/7840)/20
% 0.19e-5/200;
e=[       0         0         0         0   12.7000         0
          0         0         0   12.7000         0         0
    -5.2000   -5.2000   15.1000         0         0         0];

sde=inv(Tbar'*RT'*cde*RT*Tbar);
dd=e*sde;
dd33=2.913;epsilon33=1276;
dd31=-1.2382;epsilon33=1276;

h=diag([1/370,1/370,1/635])*e/8.85e-12;
cdd=e'*h+inv(sde);
% cdd=h'*diag([370,370,635])*h*8.85e-12+inv(sde)
epsilon_T=diag([370,370,635])+dd*e'/8.85e-12
sdd=inv(cdd)
