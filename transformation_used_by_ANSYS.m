%%%%%to transform stresses and strain from standard notation to ANSYS
%%%%%notation

clear;
clc;
syms C c11 c12 c13 c14 c15 c16 c21 c22 c23 c24 c25 c26 c31 c32 c33 c34 c35 c36
syms c41 c42 c43 c44 c45 c46 c51 c52 c53 c54 c55 c56 c61 c62 c63 c64 c65 c66
C = [c11 c12 c13 c14 c15 c16;
     c21 c22 c23 c24 c25 c26;
     c31 c32 c33 c34 c35 c36;
     c41 c42 c43 c44 c45 c46;
     c51 c52 c53 c54 c55 c56;
     c61 c62 c63 c64 c65 c66];
  
T = eye(6);
T([4,5,6],:) = T([6,4,5],:); %%%%contracted notation convention used by ANSYS

C_ANSYS = T*C/T;
disp(C_ANSYS)
