import numpy as np
#from mpmath import mpf, sqrt, power, mp
#mp.dps=100

def system01_p_l_kd__pl(p,l, kdpl):
    return (p+kdpl+l-np.lib.scimath.sqrt(-4*p*l+np.power(p+kdpl+l,2)))/2.

def system03_p_kdpp__pp(p,kdpp):
    return (4*p+kdpp-np.lib.scimath.sqrt(kdpp)*np.lib.scimath.sqrt(8*p+kdpp))/8.

# def system_HP_01_p_l_kd__pl(p,l, kdpl):
#     p=mpf(p)
#     l=mpf(l)
#     kdpl=mpf(kdpl)
#     return (p+kdpl+l-sqrt(-4*p*l+power(p+kdpl+l,2)))/2.