library(bigWig)
library(zoo)
system('wget https://github.com/guertinlab/modeling_PRO_composites/raw/main/gr.only.coordinates.Rdata')
load('gr.only.coordinates.Rdata')


system('wget https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_t0_plus_merged_normalized.bigWig')
system('wget https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_t20min_plus_merged_normalized.bigWig')
system('wget https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_t0_minus_merged_normalized.bigWig')
system('wget https://data.cyverse.org/dav-anon/iplant/home/guertin/adipo_hub/mm10/3T3_t20min_minus_merged_normalized.bigWig')


plot.composites <- function(dat, fact = 'RNA polymerase', 
                                               summit = 'TSS', class= '', num.m = -200, 
                                               num.p =90, y.low =0, y.high = 0.2, 
                                               col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  
                                               rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
                                               rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), 
                                               rgb(1/2,1/2,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), 
              rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
    pdf(paste('composite_', fact, '_signals_', summit, '_peaks', class, '.pdf', sep=''), width=8.83, 
        height=3) 
    print(xyplot(est ~ x|factor, group = time, data = dat,
                 type = 'l',
               scales=list(x=list(cex=0.8,relation = "free", axs ="i"), y =list(cex=0.8, relation="free")),
                 col = col.lines,
                 auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(strip.background=list(col="grey85"),
                                   superpose.symbol = list(pch = c(16),
                                                           col=col.lines, cex =0.5), 
                   superpose.line = list(col = col.lines, lwd=c(2), 
                       lty = c(1))),
                 cex.axis=1.0,
                 par.strip.text=list(cex=0.9, font=1, col='black'),
                 aspect=1.0,
                 between=list(y=0.5, x=0.5),
                 lwd=2,
                 ylab = list(label = paste(fact," Density", sep=''), cex =0.8),
               xlab = list(label = paste("Distance from ", summit, sep=''), cex =0.8),
               panel = function(x, y, ...) {
           panel.xyplot(x, y, ...)
           panel.abline(h = 0, lty =1, lwd = 1.0, col = 'grey85')
       }
                 ))
    dev.off()
}



composite.func <- function(gene.bed, bwPlus, bwMinus, upstream = 249,
                           downstream = 1000, step = 2, roll.avg = 10,
                           factor = 'GR', time.pt = 't0' ) {
    require(bigWig)
    bw.plus = load.bigWig(bwPlus)
    bw.minus = load.bigWig(bwMinus)
    bedTSS = fiveprime.bed(gene.bed, upstreamWindow = upstream,
                           downstreamWindow = downstream)
    tss.matrix = bed6.step.bpQuery.bigWig(bw.plus, bw.minus , bedTSS,
                                          step = step, as.matrix=TRUE, follow.strand=TRUE)
    coordin.start = (-upstream - 1)  + (step * roll.avg)/2
    coordin.end = downstream - (step * roll.avg)/2
    composite.lattice = data.frame(seq(coordin.start, coordin.end, by = step),
                                     rollmean(colMeans(tss.matrix), roll.avg),
                                     factor,
                                     time.pt,
                                     stringsAsFactors = FALSE)
    colnames(composite.lattice) = c('x', 'est', 'factor', 'time')
    composite.lattice$x = as.numeric(composite.lattice$x)
    unload.bigWig(bw.plus)
    unload.bigWig(bw.minus)
    return(composite.lattice)
}


#incorporate this into a loop

gr.0min.composite = composite.func(gene.bed = gr.only.coordinates,
               bwPlus = '3T3_t0_plus_merged_normalized.bigWig',
               bwMinus = '3T3_t0_minus_merged_normalized.bigWig',
               upstream = 499,
               downstream = 1000,
               step = 2,
               roll.avg = 20,
               factor = 'GR',
               time.pt = '0min')

gr.20min.composite = composite.func(gene.bed = gr.only.coordinates,
               bwPlus = '3T3_t20min_plus_merged_normalized.bigWig',
               bwMinus = '3T3_t20min_minus_merged_normalized.bigWig',
               upstream = 499,
               downstream = 1000,
               step = 2,
               roll.avg = 20,
               factor = 'GR',
               time.pt = '20min')

composite.lattice.gr = rbind(gr.0min.composite, gr.20min.composite)
           
plot.composites(composite.lattice.gr,col.lines = c(rgb(0,0,1), rgb(1,0,0)))
