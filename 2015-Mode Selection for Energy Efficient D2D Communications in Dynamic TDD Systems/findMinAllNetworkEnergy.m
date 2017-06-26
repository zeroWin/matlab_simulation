function f = findMinAllNetworkEnergy(x,btgt,T,W,G10,G02)
f = (exp(btgt/(W*(T-x)))-1)*(T-x)/G10 + (exp(btgt/(W*x))-1)*x/G02;
