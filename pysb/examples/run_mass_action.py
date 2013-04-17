from pysb.integrate import odesolve
from pylab import *

from pysb.examples.mass_action import model

from pysb.smoldynlib import *
from pysb.generator.smoldynlib import *





t = linspace(0, 6*60, 6*60+1)  
for a in [5, 50, 500, 5000]:
    y = odesolve(model, t, y0=[a, 5.0, 0.0])
    plot(t, y['Na']/a)

show()


# This method finds a bug!!
from pysb.integrate import odesolve
from pylab import *

from pysb.examples.mass_action import model

from pysb.smoldynlib import *
from pysb.generator.smoldynlib import *
model.parameters['Ainit'].value = 50.0
model.parameters['Binit'].value = 5.0
t = linspace(0, 60, 60+1)  # 1 hours
y = odesolve(model, t)
y['Na']

#show()  