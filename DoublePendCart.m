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
C = eye(6);
%C = [1,zeros(1,5);
%    0,0,1,zeros(1,3);
%   zeros(1,4),1,0];
%Controllability (fully controllable)
rank(ctrb(A,B));
K_c = place(A,B,[-8,-30,-2,-100,-3,-120]);

%Observability (fully observable)
rank(obsv(A,C));
L = place(A',C',10*[-8,-30,-2,-100,-3,-120])';


%% Discrete model
 dt = 0.01;
N_samples = 500;

x_lin = zeros(size(A,1),N_samples);
x0 = [0;0;pi/10;0;0;0];
x_lin(:,1) = x0;
y_lin = zeros(size(C,1),N_samples);
y_lin(:,1) = C*x0;
xhat_lin = zeros(6,N_samples);
x0hat = [0;3;pi;10;-1;5];
xhat_lin(:,1)=x0hat;

u_c_lin = zeros(1,N_samples);
u_co = zeros(1,N_samples);

%Discrete linear model
Ad = eye(6)+dt*A; Bd = dt*B;
w= 0.001; v = 1e-4;
Q = w^2*eye(size(Ad,1)); R = v^2*eye(size(C,1));
P = 1e3*eye(size(Ad,1)); P_pred = 1e10*eye(size(Ad,1));
%Discrete controllability(fully controllable)
polesd = exp([-1,-10,-1.8,-8,-2.1,-20]*dt);
Kd = place(Ad,Bd,polesd);
%Discrete observability(fully observable)

obspolesd =  exp([-1 -5 -20 -100 -30 -150]*dt);
Ld = place(Ad',C',obspolesd)';

x = zeros(6,N_samples);
x(:,1) = x0;
x_pred = zeros(size(Ad,1),N_samples);
x_pred(:,1) = zeros(size(Ad,1),1);
x_hat  = zeros(size(Ad,1),N_samples);
x_hat(:,1) = x0hat;
u_c = zeros(1,N_samples);

%%EKF
A_nl = jacobian(Ydd,y);
A_nl=[A_nl(:,1),A_nl(:,4),A_nl(:,2),A_nl(:,5),A_nl(:,3),A_nl(:,6)];
A_nl=[A_nl(1,:);A_nl(4,:);A_nl(2,:);A_nl(5,:);A_nl(3,:);A_nl(6,:)];
B_nl = jacobian(Ydd,u);
B_nl=[B_nl(1);B_nl(4);B_nl(2);B_nl(5);B_nl(3);B_nl(6)];
f_ekf(xc,xcd,alpha,alphad,theta,thetad,u) = subs(A_nl+B_nl,[a1;a2;a3;a4;a5;a6;zeta1;zeta2],[a1_c;a2_c;a3_c;a4_c;a5_c;a6_c;zeta1_c;zeta2_c]);

for i=1:N_samples-1
    %Discrete linear system
    %u_c_lin(i) = -Kd*x_hat(:,i);
    %x_lin(:,i+1) = Ad*x_lin(:,i)+Bd*u_c_lin(i)+w*randn(size(Ad,1),1);
    %y_lin(:,i+1) = C*x_lin(:,i+1)+v*randn(size(C,1),1);

    %Discrete observer
    %xhat_lin(:,i+1) = Ad*xhat_lin(:,i)+Bd*u_co(i)+Ld*(C*(x_lin(:,i)-xhat_lin(:,i)));
    %%Discrete full state feedback control
    %u_co(i+1) = -Kd*x_lin(:,i);
    
    %Kalman Filter
    
    %x_pred(:,i+1)=Ad*x_hat(:,i)+Bd*u_c_lin(i);
    %P_pred = Ad*P*Ad'+Q;
    %K_k = P_pred*C'/(C*P_pred*C'+R);
    %x_hat(:,i+1) = x_pred(:,i+1)+K_k*(y_lin(:,i+1)-C*x_pred(:,i+1));
    %P = (eye(size(Ad,1))-K_k*C)*P_pred*(eye(size(Ad,1))-K_k*C)'+K_k*R*K_k';    

    %Discrete nonlinear model
    x(:,i+1) = fd(x(:,i),-Kd*x(:,i),dt, a1_c, a2_c, a3_c, a4_c, a5_c, a6_c, zeta1_c, zeta2_c);
    %Extended Kalman Filter
    x(:,i+1) = x(:,i+1)+w*randn(size(Ad,1),1);
    y_nl(:,i+1) = C*x(:,i+1)+v*randn(size(C,1),1);
    x_pred(:,i+1) = fd(x_hat(:,i),-Kd*x(:,i),dt, a1_c, a2_c, a3_c, a4_c, a5_c, a6_c, zeta1_c, zeta2_c);
    A_ekf = eval(f_ekf(x_hat(1,i),x_hat(2,i),x_hat(3,i),x_hat(4,i),x_hat(5,i),x_hat(6,i),-Kd*x(:,i)));
    P_pred = A_ekf*P*A_ekf'+Q;
    K_ekf = P_pred*C'/(C*P_pred*C'+R);
    x_hat(:,i+1) = x_pred(:,i+1)+K_ekf*(y_nl(:,i+1)-C*x_pred(:,i+1));
    P = (eye(size(A_ekf,1))-K_ekf*C)*P_pred*(eye(size(A_ekf,1))-K_ekf*C)'+K_ekf*R*K_ekf';    
end
%% Plots

%{
subplot(321)
hold on
plot(dt*(0:N_samples-1),x_lin(1,:),dt*(0:N_samples-1),x_hat(1,:))%,dt*(0:N_samples-1),x_pred(1,:))
legend('Cart position','Cart postion(Obs.)','Cart Position(Measured)')
subplot(322)
plot(dt*(0:N_samples-1),x_lin(2,:),dt*(0:N_samples-1),x_hat(2,:))
legend('Cart velocity','Cart velocity(Pred.)')
subplot(323)
plot(dt*(0:N_samples-1),x_lin(3,:)*180/pi,dt*(0:N_samples-1),x_hat(3,:)*180/pi)%,dt*(0:N_samples-1),x_pred(1,:)*180/pi)
legend('Alpha angle','Alpha angle(Obs.)','Alpha(Measured)')
subplot(324)
plot(dt*(0:N_samples-1),x_lin(4,:)*180/pi,dt*(0:N_samples-1),x_hat(4,:)*180/pi)
legend('Alpha ang. velocity','Alpha ang. velocity(Obs.)')
subplot(325)
plot(dt*(0:N_samples-1),x_lin(5,:)*180/pi,dt*(0:N_samples-1),x_hat(5,:)*180/pi)%,dt*(0:N_samples-1),x_pred(1,:)*180/pi)
legend('Theta angle','Theta angle(Obs.)','Theta(Measured)')
subplot(326)
plot(dt*(0:N_samples-1),x_lin(6,:)*180/pi,dt*(0:N_samples-1),x_hat(6,:)*180/pi)
legend('Theta ang. velocity','Theta ang. velocity(Obs.)')
%}
%%{
figure
subplot(321)
hold on
plot(dt*(0:N_samples-1),x(1,:),dt*(0:N_samples-1),x_hat(1,:))
legend('Cart position','Cart postion(Obs.)')
subplot(322)
plot(dt*(0:N_samples-1),x(2,:),dt*(0:N_samples-1),x_hat(2,:))
legend('Cart velocity','Cart velocity(Obs.)')
subplot(323)
plot(dt*(0:N_samples-1),x(3,:)*180/pi,dt*(0:N_samples-1),x_hat(3,:)*180/pi)
legend('Alpha angle','Alpha angle(Obs.)')
subplot(324)
plot(dt*(0:N_samples-1),x(4,:)*180/pi,dt*(0:N_samples-1),x_hat(4,:)*180/pi)
legend('Alpha ang. velocity','Alpha ang. velocity(Obs.)')
subplot(325)
plot(dt*(0:N_samples-1),x(5,:)*180/pi,dt*(0:N_samples-1),x_hat(5,:)*180/pi)
legend('Theta angle','Theta angle(Obs.)')
subplot(326)
plot(dt*(0:N_samples-1),x(6,:)*180/pi,dt*(0:N_samples-1),x_hat(6,:)*180/pi)
legend('Theta ang. velocity','Theta ang. velocity(Obs.)')
figure
plot(dt*(0:N_samples-1),-Kd*x(:,:))
%}

%% Functions
function xd = fd(x,u,dt, a1_c, a2_c, a3_c, a4_c, a5_c, a6_c, zeta1_c, zeta2_c)
xd = zeros(6,1);
xd(1) = x(1)+dt*x(2);
xd(2) = x(2)+dt*((zeta1_c*sin(x(3))*(a2_c*a6_c*cos(x(3)) - a3_c*a5_c*cos(x(5))*cos(x(3) - x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) - x(6)*((a5_c*x(6)*sin(x(3) - x(5))*(a2_c*a6_c*cos(x(3)) - a3_c*a5_c*cos(x(5))*cos(x(3) - x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) + (a3_c*x(6)*sin(x(5))*(a4_c*a6_c - a5_c^2*cos(x(3) - x(5))^2))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5)))) - (u*(a4_c*a6_c - a5_c^2*cos(x(3) - x(5))^2))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) - x(4)*((a2_c*x(4)*sin(x(3))*(a4_c*a6_c - a5_c^2*cos(x(3) - x(5))^2))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) - (a5_c*x(4)*sin(x(3) - x(5))*(a3_c*a4_c*cos(x(5)) - a2_c*a5_c*cos(x(3))*cos(x(3) - x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5)))) + (zeta2_c*sin(x(5))*(a3_c*a4_c*cos(x(5)) - a2_c*a5_c*cos(x(3))*cos(x(3) - x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))));
xd(3) = x(3)+dt*x(4);
xd(4) = x(4)+dt*(x(6)*((a5_c*x(6)*sin(x(3) - x(5))*(a1_c*a6_c - a3_c^2*cos(x(5))^2))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) + (a3_c*x(6)*sin(x(5))*(a2_c*a6_c*cos(x(3)) - a3_c*a5_c*cos(x(5))*cos(x(3) - x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5)))) + x(4)*((a2_c*x(4)*sin(x(3))*(a2_c*a6_c*cos(x(3)) - a3_c*a5_c*cos(x(5))*cos(x(3) - x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) + (a5_c*x(4)*sin(x(3) - x(5))*(a1_c*a5_c*cos(x(3) - x(5)) - a2_c*a3_c*cos(x(3))*cos(x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5)))) + (u*(a2_c*a6_c*cos(x(3)) - a3_c*a5_c*cos(x(5))*cos(x(3) - x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) - (zeta1_c*sin(x(3))*(a1_c*a6_c - a3_c^2*cos(x(5))^2))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) + (zeta2_c*sin(x(5))*(a1_c*a5_c*cos(x(3) - x(5)) - a2_c*a3_c*cos(x(3))*cos(x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))));
xd(5) = x(5)+dt*x(6);
xd(6) = x(6)+dt*((u*(a3_c*a4_c*cos(x(5)) - a2_c*a5_c*cos(x(3))*cos(x(3) - x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) - x(4)*((a5_c*x(4)*sin(x(3) - x(5))*(a1_c*a4_c - a2_c^2*cos(x(3))^2))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) - (a2_c*x(4)*sin(x(3))*(a3_c*a4_c*cos(x(5)) - a2_c*a5_c*cos(x(3))*cos(x(3) - x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5)))) - x(6)*((a5_c*x(6)*sin(x(3) - x(5))*(a1_c*a5_c*cos(x(3) - x(5)) - a2_c*a3_c*cos(x(3))*cos(x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) - (a3_c*x(6)*sin(x(5))*(a3_c*a4_c*cos(x(5)) - a2_c*a5_c*cos(x(3))*cos(x(3) - x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5)))) - (zeta2_c*sin(x(5))*(a1_c*a4_c - a2_c^2*cos(x(3))^2))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5))) + (zeta1_c*sin(x(3))*(a1_c*a5_c*cos(x(3) - x(5)) - a2_c*a3_c*cos(x(3))*cos(x(5))))/(a2_c^2*a6_c*cos(x(3))^2 - a1_c*a4_c*a6_c + a3_c^2*a4_c*cos(x(5))^2 + a1_c*a5_c^2*cos(x(3) - x(5))^2 - 2*a2_c*a3_c*a5_c*cos(x(3))*cos(x(5))*cos(x(3) - x(5)))); 
end