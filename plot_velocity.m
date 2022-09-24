clearvars
close all
clc


ORA_DELLA_PARTENZA = 7.75;


lengths = [0.61, 0.23, 0.34, 0.54, 0.29, 0.31, 0.59, 0.6, 0.41, 0.2, 0.7, 0.53, 0.51];

V = readtable("C:\A_Tesi\Python\CTM-s\model\velocity.csv");
vel = table2array(V);
KONST = 360;

ti = int64(ORA_DELLA_PARTENZA * KONST);
v_plot = zeros(length(lengths), 1);
times = zeros(length(lengths)+1, 1);
times(1) = ti;
t_next=ti;
for j = 1:13
    v=vel(t_next,j);
    fanculo= int64((lengths(j)/v)*KONST);
    t_next = t_next + fanculo;
    v_plot(j) = v;
    times(j+1)=t_next; 

end
distance = 0;
x_axis = linspace(1, 24, length(v_plot));
figure(6)
hold on
grid on

for j = 1:13

  f = @(t) v_plot(j)*(t - times(j)) + distance;
  distance = lengths(j)*1000 + distance;
  fplot(f, [times(1) times(14)])
  

end 
  legend("cell 1","cell 2","cell 3","cell 4","cell 5","cell 6","cell 7","cell 8","cell 9","cell 10","cell 11","cell 12","cell 13")
for j = 1:14
  xline(times(j))

end 


