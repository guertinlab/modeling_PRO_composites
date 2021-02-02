library(lattice)

plot.pro.simulation.composites <- function(input,
                                           filename = 'dynamic_pro_model') {
    colors.ramp = colorRampPalette(c( "pink", "#ff0000", "#cd0000", "#8b0000"),
                                   bias=1, alpha = FALSE)(21)
    pdf(paste0(filename, '.pdf'), width=3.43, height = 9.43)
    print(
        xyplot(input[,2] ~ input[,3], groups = input[,1], data = input, type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"),
                   y =list(cex=0.8,alternating=FALSE)),
               col = colors.ramp,
               par.settings = list(strip.background=list(col="grey85"), 
                                   superpose.symbol = list(pch = c(16),
                                                col = colors.ramp, cex =0.7),
                           superpose.line = list(lwd=c(2.5), col = colors.ramp,
                            lty = c(1,1))),
       aspect=1.0,
#    auto.key = list(points=F, lines=T, cex=0.8),
       lwd=2,
       ylab = list(label = "Simulated PRO-seq Signal", cex =1,font=1),
       xlab = list(label = 'Distance along gene body', cex =1,font=1),
       panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.abline(h = 0, lty =1, lwd = 1.5, col = 'grey80')
       }
       )
       )
    dev.off()
}

#combining the previous three functions for a dynamic PRO-seq profile
dynamic.pro.profile <- function(input, tau = 20,min.pk = 0, max.pk =1,
                                dpk =0.001, gene.len = 1000, filename = 'dynamic_pro_model_density') {
    bp = seq(0,gene.len-1)
    colnames(input) = c("time","Pause", "Body")
    body.parameters = mapply(find.body.param,
                             bpeak = input$Body,
                             tau = tau,
                             pausepeak = input$Pause,
                             min.pk = min.pk,
                             max.pk = max.pk,
                             dpk = dpk)
    input$bodyParam = body.parameters
    df.out = as.data.frame(matrix(NA, nrow=0, ncol = 3))
    bp = seq(0,gene.len-1)
    for (i in 1:length(input$Promoter)) {#need to revisit and use apply func.
        x = get.pro.waveform(pausepeak=input$Pause[i],
                             bpeak=input$Body[i],
                             bodypeak=input$bodyParam[i],
                             bp.seq=bp,
                             tau=tau)$vec
        y = as.data.frame(cbind(i, x, bp))
        colnames(y) = colnames(df.out) 
        df.out = rbind(df.out, y)
    }
    plot.pro.simulation.composites(df.out, 
                               filename = filename) 
    return(df.out)
}

plot.changes.wrt <- function(input, filename = 'dynamic_pro_model') {
#    reformat to lattice
    lat.b = cbind(input[,c(1,3)], 'Body')
    colnames(lat.b) = c('time', 'signal', 'region')
    lat.p = cbind(input[,c(1:2)], 'Pause')
    colnames(lat.p) = colnames(lat.b)
    lattice.result = rbind(lat.p, lat.b)
  colors.ramp = colorRampPalette(c( "pink", "#ff0000", "#cd0000", "#8b0000"),
                                   bias=1, alpha = FALSE)(21)
#plot
    pdf(paste0(filename, '.pdf'), width=4.43, height = 3.43)
    print(
        xyplot(signal ~ time| region, data = lattice.result, type = c('l', 'p'),
               scales=list(x=list(cex=0.8, relation = "free"),
                   y =list(cex=0.8, relation = "free", alternating=FALSE)),
               col = 'black',
               par.settings = list(strip.background=list(col="grey85"),
                                   superpose.symbol = list(pch = c(16),
                                                           col='black', cex =0.7),
                                   superpose.line = list(col = colors.ramp, lwd=c(2.5),
                                                         lty = c(1,1,1,1,1,1,1,1,1))),
       aspect=1.0,
       between=list(y=0.5, x=0.5),
       lwd=2,
       ylab = list(label = "Density", cex =1,font=1),
       xlab = list(label = 'Time', cex =1,font=1)
       )
)
    dev.off()
}
