kinit = 1.0 #initiation rate (pooled in with the unbound polII concentration) RNA polymerases/second
krel = 0.625 #pause release rate (rate-limiting at most genes)  RNA polymerases/second
kpre = 0.375 #premature release rate (controversial how prevalent this is see David Price H2O2) RNA polymerases/second
kelong = 50 #elongation rate: 50 RNA polymerases/second RNA polymerases/second
times = seq(0,10, by = 0.5) #time

params = c(kinit=kinit, kpre=kpre, kelong=kelong, krel=krel)
initial.state = c(P = 1, B1 = 0.01, B2 = 0, B3 = 0)
 
#Solve series of differential equations
result = ode(y = initial.state,
             times = times,
             func = density.prime,
             parms = params)

result = data.frame(result)

#plot in lattice
plot.changes.wrt(result, filename = 'dynamic_pro_model_PGBdensities')


#prepare the simulated composites
dynamic.pro = dynamic.pro.profile(input = result)

#plot the steady state profile as well.
#set to zero and solve
# dP = kinit - (kpre + krel)*P
# dB = (krel)*P - kelong*B
pause.peak = kinit/(kpre + krel) 
body.peak = ((krel)*pause.peak)/kelong
steady.state.pro = dynamic.pro.profile(input = data.frame(cbind(1, pause.peak, body.peak)))

#plot composites in lattice
plot.pro.simulation.composites(dynamic.pro,
             filename = 'dynamic_pro_model_density') 

plot.pro.simulation.composites(steady.state.pro ,
             filename = 'steady_state_pro_model_density') 
