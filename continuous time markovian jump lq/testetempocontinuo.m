rng("shuffle")
%valores

T=3;
dt = 1e-4;
RS = 10;
it = T/dt;

R0 = 19.9;
a2 = 6.72e-3;
m = 22.6e-3;
g = 9.81;
L0 = 24.9e-3;
L = 0.520;

x1eq = 4.5e-3;
x2eq = 0;
x3eq = sqrt((1/L0) * m * g * (2 * a2 * (1 + (x1eq / a2))^2));

k1 = L0 * x3eq / (a2 * (1 + x1eq / a2)^2);
k2 = L0 * x3eq^2 / (a2^2 * (1 + x1eq / a2)^3);

mc=1e6;
format long;
n = 2; %número de estados da cadeia de Markov

Q = zeros(3,3,n);
R = zeros(1,n);
A = zeros(3,3,n);
B = zeros(3,1,n);

Q(:,:,1) = [1e2, 0, 0;
            0, 1, 0;
            0, 0, 1];

Q(:,:,2) = [1e2, 0, 0;
            0, 1, 0;
            0, 0, 1];

A(:,:,1) = [0, 1, 0;
            k2/m, 0, -k1/m;
            0, 0, -R0/L];

B(:,:,1) = [0;
            0;
            1/L];

A(:,:,2) = [0, 1, 0;
            k2/(2*m), 0, -k1/(2*m);
            0, 0, -R0/L];

B(:,:,2) = [0;
            0;
            1/L];

R(1,1) = 0.01;
R(1,2) = 0.01;

x = zeros(3,it/RS);
u = zeros(1,it/RS);
P = zeros(3,3,n,it);
dP = zeros(3,3,n,it); %na realidade aqui vão ser armazenados os -dPi/dt
mP = zeros(3,3,n,it);
Lambda = [ -10 , 10 ; 5 , -5 ];
p = [0 , 1;1 , 0];
G = zeros(1,3,n,it);
Ic = zeros(1,mc);

%Resolver edo (método de euler)

for k = it:-1:2
  for i = 1:n
    for j = 1:n
      mP(:,:,i,k) = mP(:,:,i,k) + Lambda(i,j)*P(:,:,j,k);
     end
    dP(:,:,i,k) = Q(:,:,i) + A(:,:,i)'*P(:,:,i,k) + P(:,:,i,k)*A(:,:,i) - P(:,:,i,k)*B(:,:,i)*inv(R(:,i))*B(:,:,i)'*P(:,:,i,k)+ mP(:,:,i,k);
    P(:,:,i,k-1) = P(:,:,i,k) + dP(:,:,i,k)*dt;
    P(:,:,i,k-1) = .5*(P(:,:,i,k-1)+P(:,:,i,k-1)');
    %[u,s,v]=svd(P(:,:,i,k-1));  s*(s>1e-8);
    %P(:,:,i,k-1) = (u*s*v');
  end
end

%simular cadeia de Markov

lim = 0;
TH = zeros(1,it/RS);

for j = 1:mc
  C = cumsum(p,n);
  y = [0];
  th = [1];



  while y(end) < T
  var = exprnd(1/-Lambda(th(end),th(end)));
  y(end+1) = y(end) + var;
  r = rand;
  th(end+1) = sum(r>C(th(end),:))+1;
  end





for k = 1:it/RS
  par = sum(k*RS>y/dt);
  TH(1,k) = th(par);
end

if TH(1,it/RS)==1
  lim = lim+1;
end

%simular as decisões u e o estado x

x(:,1) = [1e-1,0,1];

c1=0;
c2=0;
c3=0;
c4=0;


  for k=1:(it/RS)-1
    u(:,k) = -inv(R(:,TH(1,k))) * B(:,:,TH(1,k))' * P(:, :, TH(1,k),k)*x(:,k);
    dxdt = A(:,:,TH(1,k))*x(:,k) + B(:,:,TH(1,k))*u(:,k);
    x(:, k+1) = x(:,k) + dt*RS * dxdt;
    if k ==1
      c1 = x(:,k)'*Q(:,:,TH(1,k))*x(:,k) + u(:,k)'*R(:,TH(1,k))*u(:,k);
    end
    if k ==2
      c2 = x(:,k)'*Q(:,:,TH(1,k))*x(:,k) + u(:,k)'*R(:,TH(1,k))*u(:,k);
    end
    if k ==3
      c3 = x(:,k)'*Q(:,:,TH(1,k))*x(:,k) + u(:,k)'*R(:,TH(1,k))*u(:,k);
   end
    if k >=4
      c4 = x(:,k)'*Q(:,:,TH(1,k))*x(:,k) + u(:,k)'*R(:,TH(1,k))*u(:,k);
      if round((k-1)/3)== (k-1)/3
        Ic(1,j) = Ic(1,j)+ 3*(c1+3*c2+3*c3+c4)*(dt*RS)/8;
      end
      c1 = c2;
      c2 = c3;
      c3 = c4;
    end
end

if round(j/1e2)== (j/1e2)
  disp(j);
end
end

%calculando os custos

Cc = sum(Ic(1,:))/mc;
Ct = x(:,1,1)'*P(:,:,1,1)*x(:,1,1);

%Testando distribuição limite cadeia de Markov

prop1 = lim/mc;

pi0 = [1,0];

prop1lim = pi0*expm(Lambda*3);

disp(prop1);
disp(prop1lim);

%figure(1); p1 = P(1,1,1,:); plot(p1);
%figure(2); p2 = P(1,2,1,:); plot(p2);
%figure(3); p3 = P(1,1,2,:); plot(p3);
%figure(4); p4 = P(1,2,2,:); plot(p4);