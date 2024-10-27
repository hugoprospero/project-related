Q = [1e2, 0, 0;
     0, 1, 0;
     0, 0, 1];

R = 19.9;

a2 = 6.72e-3;
m = 22.6e-3;
g = 9.81;
L0 = 24.9e-3;
L = 0.520;
delT = 1e-2;
T = 10;

x = zeros(3, T/delT);
u = zeros(1, T/delT);
P = zeros(3, 3, T/delT);

x1eq = 4.5e-3;
x2eq = 0;
x3eq = sqrt((1/L0) * m * g * (2 * a2 * (1 + (x1eq / a2))^2));

k1 = L0 * x3eq / (a2 * (1 + x1eq / a2)^2);
k2 = L0 * x3eq^2 / (a2^2 * (1 + x1eq / a2)^3);

A = [0, 1, 0;
     k2/m, 0, -k1/m;
     0, 0, -R/L];

B = [0;
     0;
     1/L];

R = 0.01;

% Resolver as EDOs
H = [A, -B*inv(R)*B'; -Q, -A'];

for k = T/delT:-1:2
    XL = expm(-H*delT) * [eye(3); P(:, :, k)];
    X = XL(1:3, 1:3);
    Lambda = XL(4:6, 1:3);
    P(:, :, k-1) = Lambda * inv(X);
end

x(:, 1) = [1e-3, 0, 0];

Ic=0;

for k = 1:T/delT
    u(:,k) = -inv(R) * B' * P(:, :, k) * x(:, k);
    dxdt = A * x(:, k) + B * u(:,k);
    x(:, k+1) = x(:,k) + delT * dxdt;
    if k ==1
      c2 = x(:,k)'*Q(:,:)*x(:,k) + u(:,k)'*R(:)*u(:,k);
    endif
    if k>=2
      c1 = c2;
      c2 = x(:,k)'*Q(:,:)*x(:,k) + u(:,k)'*R*u(:,k);
      Ic = Ic + (c1+c2)*delT/2;
    end
end

Ct = x(:,1)'*P(:,:,1)*x(:,1);

