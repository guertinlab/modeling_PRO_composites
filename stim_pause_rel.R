
#initial conditions are:
#pause.peak 
#body.peak
krel.stim=0.2
result.stimulate.pause.release = ode(c(P = pause.peak , B = body.peak),
             times = times, 
             func = density.prime,
             parms = c(kinit=kinit,
                        kpre=kpre,
                        kterm=kterm,
                        krel=krel.stim,
                        lp = lp,
                        lb = lb))

result.stimulate.pause.release = data.frame(result.stimulate.pause.release)

#plot in lattice
plot.changes.wrt(result.stimulate.pause.release,
                 filename = 'dynamic_pro_model_stim_pause_release')

plot.changes.wrt(data.frame(ode(y = c(P = pause.peak , B = body.peak),
             times = seq(0,10000, by = 1),
             func = density.prime,
             parms = c(kinit=kinit,
                        kpre=kpre,
                        kterm=kterm,
                        krel=krel.stim,
                        lp = lp,
                        lb = lb))), 'dynamic_pro_model_stim_pause_PGBdensities_ss')

#plot the steady state profile as well.
pause.peak.stim = kinit/(kpre + krel.stim) 
body.peak.stim = (((lp/lb)*krel.stim)*pause.peak.stim)/kterm
steady.state.stim.pro = dynamic.pro.profile(input = data.frame(cbind(22, pause.peak.stim, body.peak.stim)))



#prepare the simulated composites
dynamic.pro.stimulate.pause.release = dynamic.pro.profile(input = result.stimulate.pause.release)

#plot composites in lattice
plot.pro.simulation.composites(dynamic.pro.stimulate.pause.release,
                               filename = 'dynamic_pro_model_density_stim_pause_release') 

#steady state composite upon stimuating pause release
plot.pro.simulation.composites(steady.state.stim.pro,
             filename = 'steady_state_pro_model_stim_pause_density') 

