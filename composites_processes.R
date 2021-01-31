
#dimensionless rates/lengths
lp = 100 #length promoter
lb = 30000 #length gene body
kinit = 0.1 #initiation rate (pooled in with the unbound polII concentration)
krel = 0.05 #pause release rate (rate-limiting at most genes) 
kpre =0.03 #premature release rate (controversial how prevalent this is see David Price H2O2)
kterm = 0.001 #termination rate relatively slow
times = seq(0,100, by = 5) #time

params = c(kinit=kinit, kpre=kpre, kterm=kterm, krel=krel, lp = lp, lb = lb)
initial.state = c(P = 0.1, B = 0.01)
 
#Solve series of differential equations
result = ode(y = initial.state,
             times = times,
             func = density.prime,
             parms = params)

result = data.frame(result)

#plot in lattice
plot.changes.wrt(result, filename = 'dynamic_pro_model_PGBdensities')

#previous plot shows that pause density reaches steady state, but 
#gene body density is still linear
plot.changes.wrt(data.frame(ode(y = initial.state,
             times = seq(0,10000, by = 1),
             func = density.prime,
             parms = params)), 'dynamic_pro_model_PGBdensities_ss')



#prepare the simulated composites
dynamic.pro = dynamic.pro.profile(input = result)

#plot the steady state profile as well.
#set to zero and solve
# dP = kinit - (kpre + krel)*P
# dB = ((lp/lb)*krel)*P - kterm*B
pause.peak = kinit/(kpre + krel) 
body.peak = (((lp/lb)*krel)*pause.peak)/kterm
steady.state.pro = dynamic.pro.profile(input = data.frame(cbind(1, pause.peak, body.peak)))

#plot composites in lattice
plot.pro.simulation.composites(dynamic.pro,
             filename = 'dynamic_pro_model_density') 

plot.pro.simulation.composites(steady.state.pro ,
             filename = 'steady_state_pro_model_density') 
