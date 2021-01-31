
library(deSolve)
library(lattice)
#import in misc. plotting and data parsing functions
source('https://raw.githubusercontent.com/guertinlab/modeling_PRO_composites/main/plotting_composites_lattice.R')

#model derived from our G&D paper (doi:10.1101/gad.328237.119):
#Declare differential equations as function
#dP/dt = kinit - (kpre +krel)p
#dB/dt = (lp/lb)*krel * p - kterm * b
#P is the first dependent variable, promoter density; dP is its derivative wrt time
#B is the second dependent variable, body density; dB is its derivative wrt time

density.prime <- function(t, initial.state, params = params)  {
    with(as.list(c(params, initial.state)), {
        dP = kinit - (kpre + krel)*P
        dB = ((lp/lb)*krel)*P - kterm*B
        res = c(dP, dB)
        list(res)
    })
}


#code adapted from our G&D paper (doi:10.1101/gad.328237.119):
#to visualize it in a composite profile form
# function for finding the gene body parameter
# bpeak: desired gene body level
# pausepeak: desired pause region level
# tau: exponential decay constant
# min.pk: minimal level for gene body peak
# max.pk: maximal level for gene body peak
# dpk: resolution for the implicit solution
find.body.param <- function(bpeak=NULL,tau=NULL,
                            pausepeak=NULL,min.pk=0,max.pk=1,
                            dpk=0.001) {
    pk = seq(min.pk,max.pk,dpk) # sequence of peak parameter values
    dat = matrix(0,length(pk),2) # data matrix
    colnames(dat) = c("body_param","body_assymp")
    for (x in 1:length(pk)) {
        bodypeak = pk[x]
        root = tau*(bodypeak+exp(1))*exp(-1)
        peak = (root/tau) * exp(-(root - tau)/tau) + bodypeak *
            (1 - exp(-root/tau))
        body = bodypeak * pausepeak / peak
        dat[x,1] = bodypeak
        dat[x,2] = body
    }
# look up the ratio of the body parameter over the assymptote
# use linear interpolation
    inter = findInterval(bpeak,dat[,2]) - 1 
    m = (dat[(inter+1),1] - dat[inter,1]) /
        (dat[(inter+1),2] - dat[inter,2])
    b = dat[inter,1] - m * dat[inter,2]
    est = m * bpeak + b
    return(est)
} 

# function to get the PRO waveform
# bpeak: desired gene body level
# pausepeak: desired pause region level
# bodypeak: body peak value to obtain a level of bpeak
# bp.seq: base pair sequence
# tau: exponential decay constant
get.pro.waveform <- function(bpeak=NULL,pausepeak=NULL,bodypeak=NULL,
                             bp.seq=NULL,tau=NULL){
    root = tau*(bodypeak+exp(1))*exp(-1)
    peak = (root/tau) * exp(-(root - tau)/tau) + bodypeak * (1 - exp(-root/tau))
    vec = sapply(bp.seq,function(x){(pausepeak/peak)*((x/tau) * exp(-(x - tau)/tau) +
                                                      bodypeak * (1 - exp(-x/tau)))})
    vec = unname(vec)
    pars = as.data.frame(rbind(cbind(max(vec),
        vec[bp.seq[length(bp.seq)]]),c(pausepeak,bpeak)))
    names(pars) = c("Peak","Assymp")
    rownames(pars) = c("simulated","requested")
    out = list(vec=vec, pars=pars)
    return(out)
}
