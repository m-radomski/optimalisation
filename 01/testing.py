from math import sqrt
from scipy.integrate import RK45
import matplotlib.pyplot as plt

def k(D):
    a = 0.98
    b = 0.63
    g = 9.81
    pp = 1

    return -a*b*D*sqrt((2*g)/pp)

def V(t, k, V0):
    x = 0.25*((k*t)**2) + k*t*sqrt(V0) + V0

    p = (-k*sqrt(V0))/(0.5*k**2)
    return x if t < p else 0

def Fin(t):
    return 0.01 * t

def Vin(t, k, V0):
    return (V0 - V(t, k, V0)) + Fin(t)

def Tin(t, k, V0):
    return ((V0 - V(t, k, V0))*90 + Fin(t)*10)/Vin(t, k, V0)

def cm2_to_m2(x):
    return x / 10**4

D = cm2_to_m2(10)

biggest_k_possible = -0.00447214

print(biggest_k_possible, k(cm2_to_m2(17)))

kiA = k(cm2_to_m2(56))
print(Vin(1, kiA, 5))

tmaxes = []
for i in range(1, 100):
    kiA = k(cm2_to_m2(i))
    kiB = k(cm2_to_m2(36.5665))

    # if(ki > biggest_k_possible):
    ts = []
    if(True):
        ode = lambda t, T: [Vin(t, kiA, 5)/(V(t, kiB, 1) + Vin(t, kiA, 5)) * (Tin(t, kiA, 5) - T)]
        last_temp = 10
        for tmax in range(2, 1001):
            sol = RK45(ode, tmax-1, [last_temp], tmax)
            while(sol.status == 'running'):
                sol.step()
            last_temp = sol.y[0]
            ts.append(sol.y[0])

    tmaxes.append(max(ts))

plt.plot(tmaxes)
plt.ylabel('some numbers')
plt.show()
