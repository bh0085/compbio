function [var,m]=my_or_modules(var1,m1,var2,m2)

var=union(var1,var2);

m=[m1,-1,m2];

[var,m]=simplify_modules(var,m);
