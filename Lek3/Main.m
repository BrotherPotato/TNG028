% I be popping bottles

s = tf("s")

G=2*(s+1)/(s^2+3*s+4)

%feedback(FG, 1) 
% is usefull if you want y feedback


%linear

p = [-2, -3]
z = -1
K = 2

G1 = zpk(z,p,K)


% ex 3.5

K2 = 0.2;
p2 = ((s^2 + s + 1)*(s + 0.2));
G2 = K2/ p2
%%
figure(3)
%a
Kp = 1
G2f = feedback(Kp * G2,1)
step(G2f)
hold on

Kp = 2
G2f = feedback(Kp * G2,1)
step(G2f)

Kp = 3
G2f = feedback(Kp * G2,1)
step(G2f)


%%
%b
clf
figure(4)
Kp = 1
G2f = feedback(Kp * G2,1)
step(G2f)
hold on

Kp = 1.5
G2f = feedback(Kp * G2,1)
step(G2f)

Kp = 2
G2f = feedback(Kp * G2,1)
step(G2f)

%%
%c
clf
figure( "Name","C")
hold on

for Kd = 0:1:3
Kp = 1
Ki = 1
%Kd = 3
T = 0.1
F = (Kp + Ki/s + Kd*s)/ (s*T+1)

G2f3 = feedback(F* G2, 1)
step(G2f3)
end

%% 2.6

z = []

Kp = 1
Pa = s^2+2*s+1
Pb = s^2 +0.4*s+1
Pc = s^2+5*s+1
Pd = s^2+s+1
Pe = s^2+2*s+4

%zpk(z, Pa, Kp)
%zpk(z, Pb, Kp)
%zpk(z, Pc, Kp)
%zpk(z, Pd, Kp)
%zpk(z, Pe, Kp)

Ga = Kp / Pa
Gb = Kp / Pb
Gc = Kp / Pc
Gd = Kp / Pd
Ge = Kp / Pe
 
figure(9)
%b
hold on
step(Ga)
pole(Ga)

step(Gb)
pole(Gb)

step(Gc)
pole(Gc)

step(Gd)
pole(Gd)

step(Ge)
plot(pole(Ge))

clf
figure("name", "AAAAA")
hold on
grid on
axis equal

plot((pole(Ga)), "xg")
plot((pole(Gb)), "xr")
plot((pole(Gc)), "xb")
plot((pole(Gd)), "xm")
plot((pole(Ge)), "xy")

figure(7)
hold on
grid on
axis equal
plot(real(pole(Ga)),imag(pole(Ga)), "xg")
plot(real(pole(Gb)),imag(pole(Gb)), "xr")
plot(real(pole(Gc)),imag(pole(Gc)), "xb")
plot(real(pole(Gd)),imag(pole(Gd)), "xm")
plot(real(pole(Ge)),imag(pole(Ge)), "xy")

%c
% further away on the LHP -> faster system
% further away on the Im axis -> more oscillations 

