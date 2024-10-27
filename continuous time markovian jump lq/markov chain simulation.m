pkg load statistics

T = 10;
format long;
y = [0];
th = [];
L = [-2,1,1;
     0,-1,1;
     10,5,-15];
P = [0,1/2,1/2;
    0,0,1;
    2/3,1/3,0];
C = cumsum(P,2);


th(end+1) = randi([1 3]);

while y(end) < T
var = exprnd(1/-L(th(end),th(end)))
  if y(end) + var<=T
    y(end+1) = y(end) + var;
    r = rand;
    th(end+1) = sum(r>C(th(end),:))+1;
  else
    y(end+1) = T
  end
end

