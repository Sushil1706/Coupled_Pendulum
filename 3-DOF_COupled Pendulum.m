% You must have a spring file with this file in the same folder

clc; 

l = 0.5; % length of the rod 
m = 2;   % mass 
a = 0.2; % first spring attached length on the rod
b = 0.4;  % % second spring attached length on the rod
k1 = 10; % first spring stiffness
k2 = 30;  % second  spring stiffness

%  theta1  = c1*a11*sin(w1*t + ph1) + c2*a21*sin(w2*t + ph2)+c3*a31*sin(w3*t + ph3)
%  theta2  = c1*a12*sin(w1*t + ph1) + c2*a22*sin(w2*t + ph2)+c3*a32*sin(w3*t + ph3)
%  theta3  = c1*a13*sin(w1*t + ph1) + c2*a23*sin(w2*t + ph2)+c3*a33*sin(w3*t + ph3)

% calcuating the c1,c2,c3,ph1,ph2,ph3 using initial condiions as aij and wi is known

%  initial conditions:-
% th1(t = 0 ) = 0
% th2(t = 0 ) = 0
% th3(t = 0 ) = x
% th1dot(t = 0 ) = 0
% th2dot(t = 0 ) = 0
% th3dot(t = 0 ) = 0
syms c1 a1 c2 a2 c3 a3 w1 w2 w3 x ; 

%  natural frequencies and relative amplitides were given in the problem 4 of assignment 4 (using that here)

w1 = 4.45;
w2 = 34.55;
w3 = 140.16;

x = 0.274;  % (initial angular displacent given to third rod in radians)

T1 =[];
T2 = [];
T3 = [];
T =[];
X1 = [];
X3 = [];
X2 = [];
Y1=[];
Y2=[];
Y3=[];
Xa1 = [];
Xa2 = [];
Xa3 = [];
Ya1 = [];
Ya2 = [];
Ya3 = [];
Xb1 = [];
Xb2 = [];
Xb3 = [];
Yb1 = [];
Yb2 = [];
Yb3 = [];


eq1 = c1*sin(a1)-1.88*c2*sin(a2)+0.0444*c3*sin(a3) == 0;
eq2 = c1*sin(a1)+0.88*c2*sin(a2)+-1.0444*c3*sin(a3) == 0;
eq3 = c1*sin(a1)+c2*sin(a2)+c3*sin(a3) == x;
eq4 = c1*w1*cos(a1)-1.88*c2*w2*cos(a2)+0.0444*c3*w3*cos(a3) == 0;
eq5 = c1*w1*cos(a1)+0.88*c2*w2*cos(a2)-1.0444*c3*w3*cos(a3) == 0;
eq6 = c1*w1*cos(a1)+c2*w2*cos(a2)+c3*w3*cos(a3) == 0;
eqns = [ eq1, eq2, eq3, eq4, eq5, eq6];
vars = [c1 c2 c3 a1 a2 a3];
[c1, c2, c3, a1, a2, a3] = vpasolve(eqns , vars);

A = 10; %time
dt = 0.5; % time interval

for i = 0:dt:A
    t = i
    
%   calculating theta1,theta2 and theta3  
th1 = 0+c1*sin(4.45*t + a1) + c2*(-1.88)*sin(34.55*t+a2)+c3*0.0444*sin(140.16*t+a3);
th2 = 0+c1*sin(4.45*t + a1) + 0.88*(c2)*sin(34.55*t+a2)+c3*(-1.0444)*sin(140.16*t+a3);
th3 = 0+c1*sin(4.45*t + a1) + (c2)*sin(34.55*t-+a2)+c3*sin(140.16*t+a3);

% calculating x and y positions of the end point and point of spring attachement of the rod
x1 = l*sin(th1);
x2 = 1 +l*sin(th2);
x3 = 2 +l*sin(th3);
y1 = -l*cos(th1);
y2 = -l*cos(th2);
y3 = -l*cos(th3);
xa1 = a*sin(th1);
xa2 = 1 + a*sin(th2);
xa3 = 2 +a*sin(th3);
ya1 = -a*cos(th1);
ya2 = -a*cos(th2);
ya3 = -a*cos(th3);
xb1 = b*sin(th1);
xb2 = 1+b*sin(th2);
xb3 = 2 + b*sin(th3);
yb1 = -b*cos(th1);
yb2 = -b*cos(th2);
yb3 = -b*cos(th3);

% storing the values
T = [T t];
T3 = [T3 th1];
T2 = [T2 th2];
T1 = [T1 th3];
X1 = [X1 x1];
X2 = [X2 x2];
X3 = [X3 x3];
Y1 = [Y1 y1];
Y2 = [Y2 y2];
Y3 = [Y3 y3];
Xa1 = [Xa1 xa1];
Xa2 = [Xa2 xa2];
Xa3 = [Xa3 xa3];
Ya1 = [Ya1 ya1];
Ya2 = [Ya2 ya2];
Ya3 = [Ya3 ya3];
Xb1 = [Xb1 xb1];
Xb2 = [Xb2 xb2];
Xb3 = [Xb3 xb3];
Yb1 = [Yb1 yb1];
Yb2 = [Yb2 yb2];
Yb3 = [Yb3 yb3];
end

B =-0.3:0.1:2.2;
c = -0.7:0.1:0.1;
ov = zeros(9);
ov1 = ones(9);
ov2 = 2*ones(9);
oh = zeros(26);

o = zeros(1,1+A/dt);
o1 = ones(1,1+A/dt);
o2 = 2*ones(1,1+A/dt);

spr = Spring(0.03,15); % Calling spring function with '0.03 m' as radius and 15 as number of coils

Vi = VideoWriter('Vibration of the 3-DoF coupled pendulum.avi'); Vi.FrameRate = 2; open(Vi); % recording the plot

for j = 1:length(T)
    set(gcf, 'Position',  [1,2, 1000,650 ]) % Setting the grahphical window size

    point1 = [Xa1(j), Ya1(j)];
    point2 = [Xa2(j), Ya2(j)];
    [a1, b1] = spr.getSpr(point1, point2); 
    plot(a1, b1,'m', 'LineWidth', 1) % plotting first spring
    hold on
    
    point3 = [Xb2(j), Yb2(j)];
    point4 = [Xb3(j), Yb3(j)];
    [a3, b4] = spr.getSpr(point3, point4);
    plot(a3, b4,'m', 'LineWidth', 1) % plotting second spring
    hold on

% plotting horizontal line on which spring is hanging 
    plot(B,oh,'k','Linewidth',1.5)
    hold on
%     plotting vertical lines to show the y-axis 
    plot(ov,c,'k', 'Linewidth',0.1)
    hold on
    plot(ov1,c,'k', 'Linewidth',0.1)
    hold on
    plot(ov2,c,'k', 'Linewidth',0.1)
    hold on

%    plotting the rods and mass at the end of the rod
    plot([o(j) Xa1(j) Xb1(j) X1(j)],[o(j) Ya1(j) Yb1(j) Y1(j)], 'r')
    hold on
    th = 0:pi/50:2*pi;    
    x1unit = 0.02 * cos(th) + X1(j);
    y1unit = 0.02 * sin(th) + Y1(j);
    x2unit = 0.02 * cos(th) + X2(j);
    y2unit = 0.02 * sin(th) + Y2(j);
    x3unit = 0.02 * cos(th) + X3(j);
    y3unit = 0.02 * sin(th) + Y3(j);
    plot([o1(j) Xa2(j) Xb2(j) X2(j)],[o(j) Ya2(j) Yb2(j) Y2(j)], 'b')
    hold on
    
    plot(x1unit, y1unit,'k','Linewidth',10);
    plot(x2unit, y2unit,'k','Linewidth',10);
    plot(x3unit, y3unit,'k','Linewidth',10);
    hold on
    plot([o2(j) Xa3(j) Xb3(j) X3(j)],[o(j) Ya3(j) Yb3(j) Y3(j)], 'g')
    hold on

    hold off
    axis([-0.5 2.5  -0.8 0.15 ])
    title('Vibration of the 3-DoF coupled pendulum','FontSize',15)
%     grid on
     xlabel('X (m)') 
     ylabel('Y (m)') 
     writeVideo(Vi, getframe(figure(1)) );
    pause(0.01)
end
close(Vi); implay('Vibration of the 3-DoF coupled pendulum.avi');

% theta vs time plot
% figure(2)
% plot(T,T1)
% xlabel('time(Sec)') 
% ylabel('theta1(rad)') 
% 
% figure(3)
% plot(T,T2)
% xlabel('time(Sec)') 
% ylabel('theta2(rad)') 
% 
% figure(4)
% plot(T,T3)
% xlabel('time(Sec)') 
% ylabel('theta3(rad)') 