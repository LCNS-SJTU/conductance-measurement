COMMENT
an synaptic current with two exponential function conductance defined from
       the book "theoretical neuroscience"
       this code is written by Songting Li and revisied by Han Hou
       
ENDCOMMENT
					       
NEURON {
	POINT_PROCESS ExpSynn
	RANGE onset, tau1, tau2, gmax, e, i ,g
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	onset=0 (ms)
	tau1=.1 (ms)	<1e-3,1e6>
	tau2=.5 (ms)
	gmax=0 	(uS)	<0,1e9>
	e=0	(mV)
}

ASSIGNED { factor (1) v (mV) i (nA)  g (uS)}

BREAKPOINT {
	if (gmax) { at_time(onset) }
    if (t<onset) {
        g = 0 
        i=0
        }
    else {
	   g = gmax * factor * (exp(-(t-onset)/tau2)-exp(-(t-onset)/tau1))
	   i = g*(v - e)
    }
}


INITIAL {
LOCAL trise
trise=(tau1*tau2)/(tau2-tau1)
factor=exp(trise/tau2*log(tau1/tau2))-exp(trise/tau1*log(tau1/tau2))
factor=1/factor
}
