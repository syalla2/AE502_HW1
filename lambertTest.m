r1 = [5644, 2830, 4170];
r2 = [-2240, 7320, 4980];
mu = 3.986004418e5;
tof=1200;
pro = 1;
ret = 0;
disp('Given in km/sec')
[v_1pro, v_2pro, ~] = LambertCurtis(r1,r2, tof, pro, mu)
[v_1ret, v_2ret, ~] = LambertCurtis(r1,r2, tof, ret, mu)