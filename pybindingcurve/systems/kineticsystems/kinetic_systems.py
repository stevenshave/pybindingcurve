import numpy as np
from scipy.integrate import solve_ivp


# 1:1 binding
def system01_p_l_kd__pl(p,l,kdpl, interval=(0, 100)):
    def ode(concs, t, kdpl):
        p, l, pl = concs
        r1 = -p*l + kdpl*pl
        dpdt = r1
        dldt = r1
        dpldt = - r1
        return [dpdt, dldt, dpldt]
    ode_result=solve_ivp(lambda t,y:ode(y,t,kdpl), interval, [p,l,0.0], rtol=1e-12, atol=1e-12).y[:,-1]
    return {'p':ode_result[0], 'l':ode_result[1],'pl':ode_result[2]}

# 1:1:1 competition
def system02_p_l_i_kdpl_kdpi__pl(p, l, i, kdpl, kdpi, interval=(0, 100)):
    def ode(concs, t, kdpl, kdpi):
        p, l, i, pl, pi  = concs
        r1 = -p*l + kdpl*pl
        r2 = -p*i + kdpi*pi
        dpdt = r1 + r2
        dldt = r1
        didt = r2
        dpldt = - r1
        dpidt = - r2
        return [dpdt, dldt, didt, dpldt, dpidt]
    ode_result=solve_ivp(lambda t,y:ode(y,t,kdpl, kdpi), interval, [p,l,i,0.0,0.0], rtol=1e-12, atol=1e-12).y[:,-1]
    return {'p':ode_result[0],'l':ode_result[1],'i':ode_result[2],'pl':ode_result[3],'pi':ode_result[4]}


# Homodimer formation
def system03_p_kdpp__pp(p,kdpp, interval=(0, 100)):
    def ode(concs, t, kdpp):
        p, pp  = concs
        r1 = -(p*p) + kdpp*pp
        dpdt = 2*r1
        dppdt = - r1
        return [dpdt, dppdt]
    ode_result=solve_ivp(lambda t,y:ode(y,t,kdpp), interval, [p,0.0], rtol=1e-12, atol=1e-12).y[:,-1]
    return {'p':ode_result[0],'pp':ode_result[1]}


# Homodlmer breaklng
def system04_p_l_kdpp_kdpl__pp(p,l,kdpp,kdpl,interval=(0, 100)):
    def ode(concs, t, kdpp, kdpl):
        p, l, pp, pl  = concs
        r_pp = -(p*p) + kdpp*pp
        r_pl = -p*l + kdpl*pl
        dpdt = 2*r_pp + r_pl
        dldt = r_pl
        dppdt = -r_pp
        dpldt = -r_pl
        return [dpdt, dldt, dppdt, dpldt]
    ode_result = solve_ivp(lambda t,y:ode(y,t,kdpp, kdpl), interval, [p, l, 0.0, 0.0], rtol=1e-12, atol=1e-12).y[:,-1]
    return {'p':ode_result[0],'l':ode_result[1],'pp':ode_result[2],'pl':ode_result[3]}

