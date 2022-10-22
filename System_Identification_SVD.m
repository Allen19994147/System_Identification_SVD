%%% System Identification
%%% Discrete system 
%%% Please read README for assumptions
%%% Allen Lee
clc
clear all
%%% Measurement data
Y = [0 1.321 -0.3612 0.336 -0.341 0.3194 -0.2994 0.2813 -0.2642...
    0.2482 -0.2332 0.219 -0.2058 0.1933 -0.1816];
y = Y(2:end);

%%% Dimension for Hankel Matrix
H_dim = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Computation %%%%%%%%%%%%%%%%%%%
H = Hankel_Matrix(y, H_dim);
[U,D,V] = svd(H);
r = rank(H);
d = sqrt(D);
P1 = [eye(r);zeros(H_dim-r,r)];
P2 = [eye(r) zeros(r,H_dim-r)];
g = d*P1;
t = P2*d;

M = U*g;
N = t*V';

M_left_inverse = P1'*inv(d)*U';
N_right_inverse = V*inv(d)*P2';
H_hat = Hankel(y(2:end),H_dim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A, B, and C matrix
B = N(:,1)
C = M(1,:)
A = M_left_inverse*H_hat*N_right_inverse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Verification part%%%%%%%%%%%%%%%%%
x0 = zeros(r,1);
y_est = [];
x = x0;

for i = 1:length(y)

    if(i==1)    % input is unit impulse
        u = 1;  % Must be consistent with experiment
    else
        u = 0;
    end
    x = A*x + B*u;
    y_est = [y_est C*x];
end

plot(y,'.',Color = 'b')
hold on
plot(y_est,'o',Color='k')
legend('Measurement','Construction')
hold off

square_error = 0;
for i=1:length(y)
    square_error = square_error + (y(i)-y_est(i))^2;
end
fprintf("The square error is %d",square_error)



