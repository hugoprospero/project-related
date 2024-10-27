#valores

T=15;
mc=1000;
format long;

Q = zeros(3,3,2);
R = zeros(1,2);
A = zeros(3,3,2);
B = zeros(3,1,2);
m = 2; #número de estados da cadeia de Markov

Q(:,:,1) = [50, 0, 0;
            0, 1, 0;
            0, 0, 1];

Q(:,:,2) = [20, 0, 0;
            0, 1, 0;
            0, 0, 1];

R(1,1) = 0.01
R(1,2) = 0.1

A(:,:,1) = [1.37, 0.02, -0.007;
            39, 1.37, -0.75;
            0, 0, 0.099];

B(:,:,1) = [0;
            -0.01;
            0.03];

A(:,:,2) = [13.7, 0.2, -0.07;
            390, 13.7, -7.5;
            0, 0, 0.099];

B(:,:,2) = [0;
            -0.1;
            3];

x = zeros(3,T,mc);
u = zeros(1,T,mc);
P = zeros(3,3,2,T);
Pmed = zeros(3,3,2,T);
p = [ .2 .8 ; .5 .5 ];
H = zeros(1,2,T);
Hi = zeros(1,2,T);
G = zeros(1,3,i,T);
c = zeros(1,mc);

#encontrar os P's

for k = T-1:-1:1,
    for i = 1:m,
      for j = 1:m,
        Pmed(:,:,i,k) = Pmed(:,:,i,k) + p(i,j)*P(:,:,j,k+1);
      end
        H(:,i,k) = R(:,i) + B(:,:,i)'*Pmed(:,:,i,k)*B(:,:,i);
        Hi(1,i,k) = inv(H(1,i,k));
        G(:,:,i,k) = -Hi(1,i,k)*B(:,:,i)'*Pmed(:,:,i,k)*A(:,:,i);
        P(:,:,i,k) = Q(:,:,i) + A(:,:,i)'*(Pmed(:,:,i,k)-Pmed(:,:,i,k)*B(:,:,i)*inv(H(:,i,k))*B(:,:,i)'*Pmed(:,:,i,k))*A(:,:,i);
    end
end

#simulando a cadeia

C = cumsum(p,2);
th = zeros(mc,T); th(:,1) = 1;
for n=1:mc,
for k=2:T,
r = rand;
th(n,k) = sum(r>C(th(n,k-1),:))+1;
end
end

#Simulando as decisões u e o estado x

for n=1:mc
x(:,1,n) = [1e-1,0,10];
end

for j=1:n,
for k=1:T-1,
    u(:,k,j) = G(:,:,th(j,k),k)*x(:,k,j);
    x(:,k+1,j) = A(:,:,th(j,k))*x(:,k,j) + B(:,:,th(j,k))*u(:,k,j);
    c(:,j) = c(:,j) + x(:,k,j)'*Q(:,:,th(j,k))*x(:,k,j) + u(:,k,j)'*R(:,th(j,k))*u(:,k,j);
end
end

Cm = sum(c(1,:))/mc;
Ct = x(:,1,1)'*P(:,:,1,1)*x(:,1,1);

