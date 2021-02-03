#initial conditions are:
#pause.peak 
#body.peak
krel.stim=krel*7
kinit.stim = kinit*4
result.stimulate.pause.release = ode(c(P = pause.peak , B1 = body.peak, B2 = body.peak, B3 = body.peak),
             times = c(0,100), 
             func = density.prime,
             parms = c(kinit=kinit.stim,
                        kpre=kpre,
                        kelong=kelong,
                        krel=krel.stim))

result.stimulate.pause.release = data.frame(result.stimulate.pause.release)
tail(result.stimulate.pause.release)

#plot in lattice
plot.changes.wrt(result.stimulate.pause.release,
                 filename = 'dynamic_pro_model_stim_pause_release')


plot.changes.wrt(data.frame(ode(y = c(P = pause.peak , B1 = body.peak, B2 = body.peak, B3 = body.peak),
             times = c(0,100),
             func = density.prime,
             parms = c(kinit=kinit.stim,
                        kpre=kpre,
                        kelong=kelong,
                        krel=krel.stim))), 'dynamic_pro_model_stim_pause_PGBdensities_ss')

#plot the steady state profile as well.
pause.peak.stim = kinit.stim/(kpre + krel.stim) 
body.peak.stim = ((krel.stim)*pause.peak.stim)/kelong
steady.state.stim.pro = dynamic.pro.profile(input = data.frame(cbind(22, pause.peak.stim, body.peak.stim)))



#prepare the simulated composites
dynamic.pro.stimulate.pause.release = dynamic.pro.profile(input = result.stimulate.pause.release,
                                                          tau = 40)

#plot composites in lattice
plot.pro.simulation.composites(dynamic.pro.stimulate.pause.release,
                               filename = 'dynamic_pro_model_density_stim_pause_release') 

#steady state composite upon stimuating pause release
plot.pro.simulation.composites(steady.state.stim.pro,
             filename = 'steady_state_pro_model_stim_pause_density') 
