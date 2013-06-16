%evaluates uniform probablilty density
function pdf=evaluateUniformDistribution(x,lower,upper)
if (x<lower || x> upper)
    pdf=0;
else
    pdf=1/(upper-lower);
end