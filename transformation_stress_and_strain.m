%%%%%In Voigt notation, the transformation matrix for the stress tensor
%%%%%or the strain tensor

clear;
clc;
syms T alpha R
syms a a11 a12 a13 a21 a22 a23 a31 a32 a33
a = [a11 a12 a13;
     a21 a22 a23;
     a31 a32 a33]; %%%%%coordinate transformation
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

disp(T)    %%stress transformation
R = diag([1 1 1 2 2 2]);  %%Reuter matrix: from tensorial shear strain to engineering strain
Tbar = R*T/R;disp(Tbar) %%strain transformation

%%%%% check inv(T)=Tbar'%%%%%%%
T0 = eval(subs(T,{a11 a12 a13 a21 a22 a23 a31 a32 a33},{cos(pi/6) sin(pi/6) 0 -sin(pi/6) cos(pi/6) 0 0 0 1}));
Tbar0 = eval(subs(Tbar,{a11 a12 a13 a21 a22 a23 a31 a32 a33},{cos(pi/6) sin(pi/6) 0 -sin(pi/6) cos(pi/6) 0 0 0 1}));

inv(T0)/Tbar0'