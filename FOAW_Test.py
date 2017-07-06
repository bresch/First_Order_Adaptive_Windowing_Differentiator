from FOAWDifferentiator import FOAWDifferentiator as FOAW
import numpy as np
import pylab as p

dt = 1 
diff = FOAW(dt,0.35)
x = [0,1,2,3,4,5,6,7,7,7,7,10,20,25,40,10,4,5,4,5,4,3,2,1,40,40,41,41,42,42,43,43,39,38,40,45]
t = range(0,len(x))

derivative = [0]*len(x)
raw_derivative = [0]*len(x)

for i in range(0,len(x)):
    derivative[i] = diff.apply(x[i]) # FOAW Derivative
    print(diff.get_last_window_size())
    if i > 1:
        raw_derivative[i] = (x[i]-x[i-1])/dt # FDM Derivative


p.step(t,x, label = 'Input signal')
p.step(t,derivative, label = 'FOAW derivative')
p.step(t,raw_derivative, label = 'FDM derivative')
p.legend()
p.show()
