from math import sqrt
import matplotlib.pyplot as plt

class Consts:
    def __init__(self, DA, DB, FRate, TA, TF):
        self.DA = DA
        self.DB = DB
        self.FRate = FRate
        self.TA = TA
        self.TF = TF

def cm2_to_m2(x):
    return x / 10**4

def k(D):
    a = 0.98
    b = 0.63
    g = 9.81
    pp = 1

    return -a*b*D*sqrt((2*g)/pp)

def dVdt(V, t, kk):
    if V > 0:
        return kk * sqrt(V)
    else:
        return 0

def Fin(t, c):
    return 0.01 * t

def fVA(VA, VB, TB, t, c):
    return dVdt(VA, t, k(c.DA))

def Fin(VA, VB, TB, t, c):
    return abs(fVA(VA, VB, TB, t, c)) + c.FRate

def Tin(VA, VB, TB, t, c):
    return ((c.TA*abs(fVA(VA, VB, TB, t, c))) + (c.TF*c.FRate)) / Fin(VA, VB, TB, t, c)

def fVB(VA, VB, TB, t, c):
    return dVdt(VB, t, k(c.DB)) + Fin(VA, VB, TB, t, c)

def fTB(VA, VB, TB, t, c):
    return Fin(VA, VB, TB, t, c)/VB * (Tin(VA, VB, TB, t, c) - TB)

tMAX = []

for DAcm2 in range(1, 101):
    c = Consts(cm2_to_m2(DAcm2), cm2_to_m2(36), 0.01, 90.0, 10.0)

    h = 1
    tfinal = 1000

    VA = [0.0] * 1000
    VB = [0.0] * 1000
    TB = [0.0] * 1000
    t = [0.0] * 1000

    VA[0] = 5.0
    VB[0] = 1.0
    TB[0] = 10.0
    t[0] = 0.0

    for i in range(tfinal - 1):
        t[i+1] = t[i] + h

        k1VA = fVA(VA[i], VB[i],TB[i],t[i], c)
        k1VB = fVB(VA[i], VB[i],TB[i],t[i], c)
        k1TB = fTB(VA[i], VB[i],TB[i],t[i], c)

        k2VA = fVA(VA[i]+(h/2)*k1VA, VB[i]+(h/2)*k1VB, TB[i]+(h/2)*k1TB, t[i]+h/2, c)
        k2VB = fVB(VA[i]+(h/2)*k1VA, VB[i]+(h/2)*k1VB, TB[i]+(h/2)*k1TB, t[i]+h/2, c)
        k2TB = fTB(VA[i]+(h/2)*k1VA, VB[i]+(h/2)*k1VB, TB[i]+(h/2)*k1TB, t[i]+h/2, c)

        k3VA = fVA(VA[i]+(h/2)*k2VA, VB[i]+(h/2)*k2VB, TB[i]+(h/2)*k2TB, t[i]+h/2, c)
        k3VB = fVB(VA[i]+(h/2)*k2VA, VB[i]+(h/2)*k2VB, TB[i]+(h/2)*k2TB, t[i]+h/2, c)
        k3TB = fTB(VA[i]+(h/2)*k2VA, VB[i]+(h/2)*k2VB, TB[i]+(h/2)*k2TB, t[i]+h/2, c)

        k4VA = fVA(VA[i]+(h  )*k3VA, VB[i]+(h  )*k3VB, TB[i]+(h  )*k3TB, t[i]+h, c)
        k4VB = fVB(VA[i]+(h  )*k3VA, VB[i]+(h  )*k3VB, TB[i]+(h  )*k3TB, t[i]+h, c)
        k4TB = fTB(VA[i]+(h  )*k3VA, VB[i]+(h  )*k3VB, TB[i]+(h  )*k3TB, t[i]+h, c)

        VA[i+1] = VA[i] + (h/6)*(k1VA + 2*k2VA + 2*k3VA + k4VA)
        VB[i+1] = VB[i] + (h/6)*(k1VB + 2*k2VB + 2*k3VB + k4VB)
        TB[i+1] = TB[i] + (h/6)*(k1TB + 2*k2TB + 2*k3TB + k4TB)

    tMAX.append(abs(max(TB) - 50))

    """
    fig, ax1 = plt.subplots()
    ax1.set_title(f"DA = {DAcm2}cm2")

    ax1.set_xlabel('Time')
    ax1.set_ylabel('Volume')
    ax1.plot(VB, label = "Volume B")
    ax1.plot(VA, label = "Volume A")
    ax1.legend(bbox_to_anchor=(0,0), loc="lower left",  bbox_transform=fig.transFigure)
    ax1.tick_params(axis='y')

    ax2 = ax1.twinx()

    ax2.set_ylabel('Temp')
    ax2.plot(TB, 'r--', label="Temperature B")
    ax2.legend(bbox_to_anchor=(1,0), loc="lower right",  bbox_transform=fig.transFigure)
    ax2.tick_params(axis='y')

    fig.tight_layout()
    plt.savefig(f'fig{DAcm2}.png')

    print(f"Done with {DAcm2}")
    """

plt.plot(tMAX)
plt.show()
