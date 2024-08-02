import numpy as np
from scipy.integrate import solve_ivp

# Define parameters
prm = [1.0, 0.5, 0.3, 0.7, 0.2, 0.1]

# Initial conditions: [s, e, c_1, c_2, p]
y0 = [1.0, 0.5, 0.0, 0.0, 0.0]

# Time span
t_span = (0, 10)
t_eval = np.linspace(0, 10, 100)

# Solve the system
def dydt_p(t, y, prm):
    # parameters
    k1 = prm[0]
    k_1 = prm[1]
    k2 = prm[2]
    k3 = prm[3]
    k_3 = prm[4]
    k4 = prm[5]

    # variables
    s = y[0]
    e = y[1]
    c_1 = y[2]
    c_2 = y[3]
    p = y[4]

    j1 = k1 * s * e
    j_1 = k_1 * c_1
    j2 = k2 * c_1
    j3 = k3 * s * c_1
    j_3 = k_3 * c_2
    j4 = k4 * c_2

    dydt = np.zeros(5)
    dydt[0] = -j1 + j_1 - j3 + j_3
    dydt[1] = -j1 + j_1 + j2
    dydt[2] = j1 - j_1 - j2 - j3 + j_3 + j4
    dydt[3] = j3 - j_3 - j4
    dydt[4] = j2 + j4

    print(f"t: {t}, dydt: {dydt}")

    return dydt

sol = solve_ivp(lambda t, y: dydt_p(t, y, prm), t_span, y0, t_eval=t_eval)