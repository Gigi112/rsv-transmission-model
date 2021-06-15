syms x y
eqn1 = x/(x+y) == 0.054; % Mean
eqn2 = (x-1/3)/(x + y -2/3) == 0.027; % Median
eqn3 = (x-1)/(x + y -2) == 0.016; %Mode

sol = solve([eqn1, eqn2], [x, y]);



sol = solve([eqn1, eqn3], [x, y]);

sol = solve([eqn2, eqn3], [x, y]);
sol.x % 719/375
sol.y % 7177/125