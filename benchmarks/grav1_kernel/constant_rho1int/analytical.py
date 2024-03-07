import numpy as np

G = 6.6723e-11
pi = np.pi
def analytical(k1,k2,k3, Rcmb, x,y,z ):
    Rcmb3 = Rcmb ** 3
    r = (x ** 2 + y ** 2 + z ** 2) ** 0.5
    r3 = r**3
    return  (4/3) * pi * G * (k1*x + k2*y + k3*z) * ((Rcmb3/r3)-1)

def analytical_alpha(ki, Rcmb,R, x,y,z ):
    Rcmb3 = Rcmb ** 3
    r = (x ** 2 + y ** 2 + z ** 2) ** 0.5
    r2 = r**2

    return -(2/3)*np.pi*G * ki *  (3*(R**2) - r2 - (2*Rcmb3/r) )



def analytical_pot_homo_fullsphere(dimrho, R, x,y,z):
    r2 = x**2 + y**2 + z**2
    R2 = R**2

    return (2/3)*np.pi * G * dimrho * (3*R2 - r2)