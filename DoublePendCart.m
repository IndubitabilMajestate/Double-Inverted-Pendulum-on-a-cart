clearvars
syms xc xcd xcdd alpha alphad alphadd theta thetad thetadd a1 a2 a3 a4 a5 a6 zeta1 zeta2 u

% General equation of motion of the system
% H(y)*yd = -M(y)*y-F*(y)+L*u
H(alpha,theta) = ...
    [           a1       a2*cos(alpha)       a3*cos(theta);
     a2*cos(alpha)                  a4 a5*cos(alpha-theta);
     a3*cos(theta) a5*cos(alpha-theta)                 a6];
M(alpha,theta,alphad,thetad) = ...
    [0       -a2*sin(alpha)*alphad      -a3*sin(theta)*thetad;
     0                           0 a5*sin(alpha-theta)*thetad;
     0 -a5*sin(alpha-theta)*alphad                          0];
F(alpha,theta) = ...
    [               0; 
    -zeta1*sin(alpha); 
    -zeta2*sin(theta)];
L = [1; 
     0; 
     0];
y = [  xc;
    alpha;
    theta;
      xcd; 
   alphad; 
   thetad];
yd = [  xcd;
     alphad;
     thetad;
       xcdd;
    alphadd;
    thetadd];

eq = yd == [eye(3),zeros(3);zeros(3),-inv(H(alpha,theta))*M(alpha,theta,alphad,thetad)]*y+[zeros(3,1);inv(H(alpha,theta))*L]*u+[zeros(3,1);-inv(H(alpha,theta))*F(alpha,theta)];

S = solve(eq,[xcd,alphad,thetad,xcdd,alphadd,thetadd]);
%%
Ydd = [zeros(3),eye(3);zeros(3),-inv(H(alpha,theta))*M(alpha,theta,alphad,thetad)]*y+[zeros(3,1);inv(H(alpha,theta))*L]*u+[zeros(3,1);-inv(H(alpha,theta))*F(alpha,theta)];
x0 = [0,0,0,0,0,0]';
mc_c = 1; m1_c = 1; m2_c = 0.5;
l1_c = 0.5; l2_c = 0.25; g_c = 9.81;
a1_c = mc_c+m1_c+m2_c; a2_c = l1_c*(m1_c/2+m2_c);a3_c = m2_c*l2_c/2; a4_c = l1_c^2*(m1_c/3+m2_c); a5_c = m2_c*l1_c*l2_c/2; a6_c = m2_c*l2_c^2/3;
zeta1_c = l1_c*g_c*(m1_c/2+m2_c); zeta2_c = m2_c*l2_c*g_c/2;

% Linearization around the UP-UP position
A =jacobian(Ydd,y);
B =jacobian(Ydd,u);
A  = subs(A,[y;u;a1;a2;a3;a4;a5;a6;zeta1;zeta2],[x0;0;a1_c;a2_c;a3_c;a4_c;a5_c;a6_c;zeta1_c;zeta2_c]);
B  = subs(B,[y;u;a1;a2;a3;a4;a5;a6;zeta1;zeta2],[x0;0;a1_c;a2_c;a3_c;a4_c;a5_c;a6_c;zeta1_c;zeta2_c]);

for i = 1:6
    for j = 1:6
        if hasSymType(A(i,j),'eq')
            A(i,j) = rhs(A(i,j));
        end
    end
    if hasSymType(B(i),'eq')
        B(i) = rhs(B(i));
    end
end
A = eval(A);
B = eval(B);
% rearrange matrix A,B
A=[A(:,1),A(:,4),A(:,2),A(:,5),A(:,3),A(:,6)];
A=[A(1,:);A(4,:);A(2,:);A(5,:);A(3,:);A(6,:)];
B=[B(1);B(4);B(2);B(5);B(3);B(6)];
C = [1,zeros(1,5);
    0,0,1,zeros(1,3);
    zeros(1,4),1,0];
%Controllability (fully controllable)
rank(ctrb(A,B));
K_c = place(A,B,[-8,-30,-2,-100,-3,-120]);

%Observability (fully observable)
rank(obsv(A,C));
L = place(A',C',10*[-8,-30,-2,-100,-3,-120])';


%% Discrete model
dt = 0.0025;
N_samples = 1000;

x_lin = zeros(6,N_samples);
x0 = [0;0;pi/5;0;0;0];
x_lin(:,1) = x0;
y_lin = zeros(3,N_samples);

xhat_lin = zeros(6,N_samples);
x0hat = [0;0;0;0;0;0];
xhat_lin(:,1)=x0hat;

u_c_lin = zeros(1,N_samples);
u_co = zeros(1,N_samples);

%Discrete linear model
Ad = eye(6)+dt*A; Bd = dt*B;
w= 0.00001; v = 0.0001;
Q = w*eye(6); R = v*eye(3);
P = eye(6); P_pred = 1e5*eye(6);

%Discrete controllability(fully controllable)
polesd = exp([-8,-30,-2,-100,-3,-120]*dt);
Kd = place(Ad,Bd,polesd);
%Discrete observability(fully observable)
obspolesd =  exp([-10 -50 -20 -100 -30 -150]*dt);
Ld = place(Ad',C',obspolesd)';

x = zeros(6,N_samples);
x(:,1) = x0;
x_pred = zeros(6,N_samples);
x_pred(:,1) = [0;0;0;0;0;0];
u_c = zeros(1,N_samples);
for i=1:N_samples-1
    %Discrete linear system
    x_lin(:,i+1) = Ad*x_lin(:,i)+Bd*u_c_lin(i)+w*randn(6,1);
    y_lin(:,i+1) = C*x_lin(:,i)+v*randn(3,1);
    %Discrete full state feedback control
    u_c_lin(i+1) = -Kd*xhat_lin(:,i);


    %Discrete observer
    %xhat_lin(:,i+1) = Ad*xhat_lin(:,i)+Bd*u_co(i)+Ld*(C*(x_lin(:,i)-xhat_lin(:,i)));
    %%Discrete full state feedback control
    %u_co(i+1) = -Kd*x_lin(:,i);
    
    %Kalman Filter
    K = (P_pred*C')/(C*P*C'+R);
    xhat_lin(:,i+1) = x_pred(:,i)+K*(y_lin(:,i)-C*x_pred(:,i));
    P = (eye(6)-K*C)*P_pred;
    x_pred(:,i+1) = Ad*xhat_lin(:,i+1)+Bd*u_c_lin(i);
    P_pred = Ad*P*Ad'+Q;

    %Wrap around [-pi,pi)
    %while x(3,i) <= -pi
    %    x(3,i) = x(3,i) + 2*pi;
    %end
    %while x(3,i) > pi
    %    x(3,i) = x(3,i) - 2*pi;
    %end
    %while x(5,i) <= -pi
    %    x(5,i) = x(5,i) + 2*pi;
    %end
    %while x(5,i) > pi
    %    x(5,i) = x(5,i) - 2*pi;
    %end
    %%Discrete nonlinear model
    %x(1,i+1) = x(1,i)+dt*x(2,i);
    %x(2,i+1) = x(2,i)+dt*((a2_c*a6_c*cos(x(3,i))*sin(x(3,i)) - a3_c*a5_c*sin(x(3,i))*cos(x(5,i))*cos(x(3,i) - x(5,i)))/(a6_c*a2_c^2*cos(x(3,i))^2 - 2*a2_c*a3_c*a5_c*cos(x(3,i))*cos(x(5,i))*cos(x(3,i) - x(5,i)) + a4_c*a3_c^2*cos(x(5,i))^2 + a1_c*a5_c^2*cos(x(3,i) - x(5,i))^2 - a1_c*a4_c*a6_c))*zeta1_c - (a4_c*a6_c*u_c(i) - a5_c^2*u_c(i)*cos(x(3,i) - x(5,i))^2 + a2_c*a4_c*a6_c*x(3,i)^2*sin(x(3,i)) + a3_c*a4_c*a6_c*x(5,i)^2*sin(x(5,i)) - a3_c*a4_c*zeta2_c*cos(x(5,i))*sin(x(5,i)) - a2_c*a5_c^2*x(3,i)^2*sin(x(3,i))*cos(x(3,i) - x(5,i))^2 - a3_c*a5_c^2*x(5,i)^2*sin(x(5,i))*cos(x(3,i) - x(5,i))^2 + a2_c*a5_c^2*x(3,i)^2*sin(x(3,i) - x(5,i))*cos(x(3,i))*cos(x(3,i) - x(5,i)) - a3_c*a4_c*a5_c*x(3,i)^2*sin(x(3,i) - x(5,i))*cos(x(5,i)) + a2_c*a5_c*a6_c*x(5,i)^2*sin(x(3,i) - x(5,i))*cos(x(3,i)) - a3_c*a5_c^2*x(5,i)^2*sin(x(3,i) - x(5,i))*cos(x(5,i))*cos(x(3,i) - x(5,i)) + a2_c*a5_c*zeta2_c*cos(x(3,i))*sin(x(5,i))*cos(x(3,i) - x(5,i)))/(a6_c*a2_c^2*cos(x(3,i))^2 - 2*a2_c*a3_c*a5_c*cos(x(3,i))*cos(x(5,i))*cos(x(3,i) - x(5,i)) + a4_c*a3_c^2*cos(x(5,i))^2 + a1_c*a5_c^2*cos(x(3,i) - x(5,i))^2 - a1_c*a4_c*a6_c);
    %x(3,i+1) = x(3,i)+dt*x(4,i);
    %x(4,i+1) = x(4,i)+dt*((a3_c^2*sin(x(3,i))*cos(x(5,i))^2 - a1_c*a6_c*sin(x(3,i)))/(a6_c*a2_c^2*cos(x(3,i))^2 - 2*a2_c*a3_c*a5_c*cos(x(3,i))*cos(x(5,i))*cos(x(3,i) - x(5,i)) + a4_c*a3_c^2*cos(x(5,i))^2 + a1_c*a5_c^2*cos(x(3,i) - x(5,i))^2 - a1_c*a4_c*a6_c))*zeta1_c + (a2_c*a6_c*u_c(i)*cos(x(3,i)) + a1_c*a5_c^2*x(3,i)^2*sin(x(3,i) - x(5,i))*cos(x(3,i) - x(5,i)) + a1_c*a5_c*a6_c*x(5,i)^2*sin(x(3,i) - x(5,i)) - a3_c*a5_c*u_c(i)*cos(x(5,i))*cos(x(3,i) - x(5,i)) + a1_c*a5_c*zeta2_c*sin(x(5,i))*cos(x(3,i) - x(5,i)) - a3_c^2*a5_c*x(5,i)^2*sin(x(3,i) - x(5,i))*cos(x(5,i))^2 + a2_c^2*a6_c*x(3,i)^2*cos(x(3,i))*sin(x(3,i)) + a2_c*a3_c*a6_c*x(5,i)^2*cos(x(3,i))*sin(x(5,i)) - a3_c^2*a5_c*x(5,i)^2*cos(x(5,i))*sin(x(5,i))*cos(x(3,i) - x(5,i)) - a2_c*a3_c*zeta2_c*cos(x(3,i))*cos(x(5,i))*sin(x(5,i)) - a2_c*a3_c*a5_c*x(3,i)^2*sin(x(3,i) - x(5,i))*cos(x(3,i))*cos(x(5,i)) - a2_c*a3_c*a5_c*x(3,i)^2*sin(x(3,i))*cos(x(5,i))*cos(x(3,i) - x(5,i)))/(a6_c*a2_c^2*cos(x(3,i))^2 - 2*a2_c*a3_c*a5_c*cos(x(3,i))*cos(x(5,i))*cos(x(3,i) - x(5,i)) + a4_c*a3_c^2*cos(x(5,i))^2 + a1_c*a5_c^2*cos(x(3,i) - x(5,i))^2 - a1_c*a4_c*a6_c);
    %x(5,i+1) = x(5,i)+dt*x(6,i);
    %x(6,i+1) = x(6,i)+dt*((a1_c*a5_c*sin(x(3,i))*cos(x(3,i) - x(5,i)) - a2_c*a3_c*cos(x(3,i))*sin(x(3,i))*cos(x(5,i)))/(a6_c*a2_c^2*cos(x(3,i))^2 - 2*a2_c*a3_c*a5_c*cos(x(3,i))*cos(x(5,i))*cos(x(3,i) - x(5,i)) + a4_c*a3_c^2*cos(x(5,i))^2 + a1_c*a5_c^2*cos(x(3,i) - x(5,i))^2 - a1_c*a4_c*a6_c))*zeta1_c + (a3_c*a4_c*u_c(i)*cos(x(5,i)) - a1_c*a4_c*zeta2_c*sin(x(5,i)) + a2_c^2*zeta2_c*cos(x(3,i))^2*sin(x(5,i)) + a3_c^2*a4_c*x(5,i)^2*cos(x(5,i))*sin(x(5,i)) - a1_c*a4_c*a5_c*x(3,i)^2*sin(x(3,i) - x(5,i)) - a2_c*a5_c*u_c(i)*cos(x(3,i))*cos(x(3,i) - x(5,i)) - a1_c*a5_c^2*x(5,i)^2*sin(x(3,i) - x(5,i))*cos(x(3,i) - x(5,i)) + a2_c^2*a5_c*x(3,i)^2*sin(x(3,i) - x(5,i))*cos(x(3,i))^2 - a2_c^2*a5_c*x(3,i)^2*cos(x(3,i))*sin(x(3,i))*cos(x(3,i) - x(5,i)) + a2_c*a3_c*a4_c*x(3,i)^2*sin(x(3,i))*cos(x(5,i)) + a2_c*a3_c*a5_c*x(5,i)^2*sin(x(3,i) - x(5,i))*cos(x(3,i))*cos(x(5,i)) - a2_c*a3_c*a5_c*x(5,i)^2*cos(x(3,i))*sin(x(5,i))*cos(x(3,i) - x(5,i)))/(a6_c*a2_c^2*cos(x(3,i))^2 - 2*a2_c*a3_c*a5_c*cos(x(3,i))*cos(x(5,i))*cos(x(3,i) - x(5,i)) + a4_c*a3_c^2*cos(x(5,i))^2 + a1_c*a5_c^2*cos(x(3,i) - x(5,i))^2 - a1_c*a4_c*a6_c);
    %Discrete nonlinear full state feedback control
    %u_c(i+1) = -Kd*x(:,i);
end
%%
subplot(321)
hold on
plot(dt*(0:N_samples-1),x_lin(1,:),dt*(0:N_samples-1),xhat_lin(1,:))
legend('Cart position','Cart postion(Obs.)')
subplot(322)
plot(dt*(0:N_samples-1),x_lin(2,:),dt*(0:N_samples-1),xhat_lin(2,:),dt*(0:N_samples-1),y_lin(1,:))
legend('Cart velocity','Cart velocity(Pred.)','Cart Velocity(Measured)')
subplot(323)
plot(dt*(0:N_samples-1),x_lin(3,:)*180/pi,dt*(0:N_samples-1),xhat_lin(3,:)*180/pi)
legend('Alpha angle','Alpha angle(Obs.)')
subplot(324)
plot(dt*(0:N_samples-1),x_lin(4,:)*180/pi,dt*(0:N_samples-1),xhat_lin(4,:)*180/pi,dt*(0:N_samples-1),y_lin(1,:)*180/pi)
legend('Alpha ang. velocity','Alpha ang. velocity(Obs.)','Alpha ang. velocity(Measured)')
subplot(325)
plot(dt*(0:N_samples-1),x_lin(5,:)*180/pi,dt*(0:N_samples-1),xhat_lin(5,:)*180/pi)
legend('Theta angle','Theta angle(Obs.)')
subplot(326)
plot(dt*(0:N_samples-1),x_lin(6,:)*180/pi,dt*(0:N_samples-1),xhat_lin(6,:)*180/pi,dt*(0:N_samples-1),y_lin(1,:)*180/pi)
legend('Theta ang. velocity','Theta ang. velocity(Obs.)','Theta ang. velocity(Measured)')
figure
plot(dt*(0:N_samples-1),y_lin)
%%
%figure
%subplot(321)
%hold on
%plot(dt*(0:N_samples-1),x(1,:),dt*(0:N_samples-1),xhat_lin(1,:))
%legend('Cart position','Cart postion(Obs.)')
%subplot(322)
%plot(dt*(0:N_samples-1),x(2,:),dt*(0:N_samples-1),xhat_lin(2,:))
%legend('Cart velocity','Cart velocity(Obs.)')
%subplot(323)
%plot(dt*(0:N_samples-1),x(3,:)*180/pi,dt*(0:N_samples-1),xhat_lin(3,:)*180/pi)
%legend('Alpha angle','Alpha angle(Obs.)')
%subplot(324)
%plot(dt*(0:N_samples-1),x(4,:)*180/pi,dt*(0:N_samples-1),xhat_lin(4,:)*180/pi)
%legend('Alpha ang. velocity','Alpha ang. velocity(Obs.)')
%subplot(325)
%plot(dt*(0:N_samples-1),x(5,:)*180/pi,dt*(0:N_samples-1),xhat_lin(5,:)*180/pi)
%legend('Theta angle','Theta angle(Obs.)')
%subplot(326)
%plot(dt*(0:N_samples-1),x(6,:)*180/pi,dt*(0:N_samples-1),xhat_lin(6,:)*180/pi)
%legend('Theta ang. velocity','Theta ang. velocity(Obs.)')
%figure
%plot(dt*(0:N_samples-1),u_c)