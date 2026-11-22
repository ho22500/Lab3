%% Lab 3 Step 1 Starter Code

clear all;
close all;
clc;
% load csv data. Put data into column vectors.
% convert from degrees to radians
M = readmatrix("mocap_data.csv");

% Number of rows in M = Number of samples from motion capture system
N = size(M,1);

t = M(:,1);
x = M(:,2);
y = M(:,3);
z = M(:,4);

p = [x y z]; % position vector (in m)

roll = deg2rad(M(:,5));
pitch = deg2rad(M(:,6));
yaw = deg2rad(M(:,7));
%% 
% *Problem 1a*
% 
% For this problem, you need to use |subplot.|

% plot x,y,z data in subplots
figure(1)
subplot(3,1,1)
plot(M(:,1), M(:,2))
title("Position vs. Time")
ylabel("X (cm)")
subplot(3,1,2)
plot(M(:,1), M(:,3))
ylabel("Y (cm)")
subplot(3,1,3)
plot(M(:,1), M(:,4))
ylabel("Z (cm)")
xlabel("time (s)")
%% 
% *Problem 1b*
% 
% We approximate the velocities using the definition of a derivative, and replacing 
% $h$ with $\Delta t$.
% 
% $$v(t) = \lim_{h \rightarrow 0^+} \frac{p(t+h)-p(t)}{h} \Rightarrowv(t_i) 
% \approx  \frac{p(t_{i+1})-p(t_i)}{t_{i+1}-t_i} $$

% Calculate velocities numerically
v = zeros(N,3);
for i=1:N-1
    v(i,:) = (M(i+1,2:4)-M(i,2:4))/0.01;
end
v(end,:) = v(end-1,:);

% plot the figure
figure(2)
subplot(3,1,1)
plot(M(:,1), v(:,1))
title("Velocity vs. Time")
ylabel("V_{x} (m/s)")
subplot(3,1,2)
plot(M(:,1), v(:,2))
ylabel("V_{y} (m/s)")
subplot(3,1,3)
plot(M(:,1), v(:,3))
ylabel("V_{z} (m/s)")
xlabel("time (s)")
%% 
% *Problem 1c*
% 
% First we need to convert roll-pitch-yaw to rotation matrices. You can either 
% use change of basis (MLS p. 31-33) or use the |eul2rotm| function from the Robotics 
% Toolbox to do the job. (This part is done for you.)

% convert rpy to rotation matrices by hand (MLS p. 31-33)
R = zeros(3,3,N);

% using inline functions to make our code easier to read
Rx = @(phi) [1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
Ry = @(beta) [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
Rz = @(alpha) [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

for i=1:N
    R(:,:,i) = Rz(yaw(i))*Ry(pitch(i))*Rx(roll(i));
end

% convert rpy to rotation matrices using Robotics Toolbox
% this code works identically to above!
% eul = [yaw pitch roll];
% Rotm = eul2rotm(eul,'ZYX');
%% 
% Now implement ECE569_MatrixLog3, ECE569_so3ToVec, ECE569_VecToso3, ECE569_AxisAng3 
% in their respective .m files. When you have written them, run the command "runtests" 
% in the Command Window.
% 
% Once you have done that, take a matrix logarithm to obtain $$\widehat{\omega} 
% \theta$ (the exponential coordinate of $R=e^{\widehat{\omega}\theta}$)$. Then, 
% perform manipulations to obtain $\hat \omega$ and $\theta$ where $\omega$ is 
% a unit vector representing the axis of rotation.

% get (w,theta) from e^(w_hat_theta)
w = zeros(N,3);
theta = zeros(N,1);
for i=1:N
    Ri = R(:,:,i);
    % Calculate the matrix logarithm as we learned in class
    % See MR textbook eqn. 3.58
    w_theta = ECE569_MatrixLog3(Ri);
    % you will need to assign values to w(i,:) and theta(i)
    [w_temp, theta(i)] = ECE569_AxisAng3(w_theta);
    w(i,:) = ECE569_so3ToVec(w_temp);
end
%% 
% Finally, plot the desired signals with appropriate labels.

% subplots
figure(3)
subplot(4,1,1)
plot(M(:,1), w(:,1))
title("Orientation vs. Time")
ylabel("\omega_{1}")
subplot(4,1,2)
plot(M(:,1), w(:,2))
ylabel("\omega_{2}")
subplot(4,1,3)
plot(M(:,1), w(:,3))
ylabel("\omega_{3}")
subplot(4,1,4)
plot(M(:,1), theta(:))
ylabel("\theta (rad)")
xlabel("time (s)")
%% 
% *Problem 1d*
% 
% Take a matrix logarithm of $R(t_i)^\top R(t_{i+1})$ and then do manipulations 
% to obtain $\omega^b(t_i)$, which we defined as $\omega^b(t_i)=\frac{\omega(t_i)\theta(t_i)}{t_{i+1}-t_i}$.
w_hat_theta = zeros(3,3,N);
w = zeros(N,3);
theta = zeros(N);
w_b = zeros(N,3);
for i=1:N-1
    w_theta = ECE569_MatrixLog3(R(:,:,i)' * R(:,:,i+1));
    [w_temp, theta(i)] = ECE569_AxisAng3(w_theta);
    w(i,:) = ECE569_so3ToVec(w_temp);
    w_b(i,:) = w(i,:)*theta(i)/0.01;
end
w_b(end,:) = w_b(end-1,:);
% Calculate angular velocity
% make subplot
figure(4)
subplot(3,1,1)
plot(M(:,1), w_b(:,1))
title("Angular Velocity vs. Time")
ylabel("w^{b}_{1} (rad/s)")
subplot(3,1,2)
plot(M(:,1), w_b(:,2))
ylabel("w^{b}_{2} (rad/s)")
subplot(3,1,3)
plot(M(:,1), w_b(:,3))
ylabel("w^{b}_{3} (rad/s)")
xlabel("time (s)")
%% 
% *Problem 1e*
% 
% Plot the function $d(R(t),R_d)^2 = 2Tr(I-R^\top R_d)$ 

% calculate d^2
Rd = R(:,:,end);
d_square = zeros(N);
for i=1:N
    d_square(i) = 2*trace(eye(3)-R(:,:,i)'*Rd);
end
% plot the error from current rotation to desired rotation Rd = R(t_N)
figure(5)
plot(M(:,1), d_square)
title("Squared Distance to Target Orientation")
ylabel("d^{2}")
xlabel("time (s)")