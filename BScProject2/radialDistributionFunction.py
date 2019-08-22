import numpy as np
import matplotlib.pyplot as plt

with open('C:\\Users\\Guy\\CLionProjects\\BScProject2\\cmake-build-debug\\radialdistfunc.txt') as f:
    lines = f.readlines()
    i = lines[0::2]
    gr1 = lines[1::2]
        
with open('C:\\Users\\Guy\\CLionProjects\\BScProject2\\cmake-build-debug\\radial2.txt') as g:
    lines1 = g.readlines()
    i2 = lines1[0::2]
    gr2 = lines1[1::2]

with open('C:\\Users\\Guy\\CLionProjects\\BScProject2\\cmake-build-debug\\radial3.txt') as h:
    lines2 = h.readlines()
    i3 = lines2[0::2]
    gr3 = lines2[1::2]
    
with open('C:\\Users\\Guy\\CLionProjects\\BScProject2\\cmake-build-debug\\radial4.txt') as j:
    lines3 = j.readlines()
    i4 = lines3[0::2]
    gr4 = lines3[1::2]
    
with open('C:\\Users\\Guy\\CLionProjects\\BScProject2\\cmake-build-debug\\radial5.txt') as k:
    lines4 = k.readlines()
    i5 = lines4[0::2]
    gr5 = lines4[1::2]  
  

fig = plt.figure(1)
ax1 = fig.add_subplot(111)
ax1.set_title("Radial Distribution Function")    
ax1.set_xlabel(r'$r$')
ax1.set_ylabel(r'g(r)')


ax1.plot(i,gr1, c='r',label = r'$\rho$ = 0.2')
ax1.plot(i2,gr2, c='b', label = r'$\rho$ = 0.5')
ax1.plot(i3,gr3, c='orange', label = r'$\rho$ = 0.7')
ax1.plot(i5,gr5, c='purple', label = r'$\rho$ = 1.2')

ax1.legend(loc='lower right', frameon=True)


plt.show()