from liblibra_core import *
import libra_py
import libra_py.models.Esch_Levine as Esch_Levine
import libra_py.models.Holstein as Holstein
import libra_py.models.GLVC as GLVC
import libra_py.units as units


def custom_GLVC_set(E1, E2, V01, V02, V12):
    """
    3-state spin-boson model: the bath parameters are from:

    Bondarenko, A. S.; Tempelaar, R. Overcoming Positivity Violations for Density Matrices in Surface Hopping. 
    The Journal of Chemical Physics 2023, 158 (5), 054117. https://doi.org/10.1063/5.0135456.

    but the number of discretized bath modes is reduced to 25

    and the system parameters are adjusted

    Args:
        E1, E2, V01, V02, V12: energies and couplings in terms of kT


    Returns:
        dictionary: params, will contain the parameters:

        * **nstates** (int): the number of electronic states

        * **num_osc** (int): the number of oscillators per electronic state

        * **spectral_density** (int): type of spectral density calculation - Debye

        * **Omega** (double): characteristic frequency of bath [units: Ha]

        * **lambda** (double): bath reorganization energy [units: Ha]

        * **omega** (list of doubles): frequencies of the discretized spectral density [ units: Ha ]

        * **coupl** (list of doubles): electron-nuclear couplings, same for all states, different for each oscillator [ units: Ha ]

        * **mass** (list of doubles): masses of nuclear DOFs [ units: a.u. ]

        * **coupling_scaling** (list of `nstates` doubles): the scaling applied to the overall linear coupling terms for each state [ unitS: N/A ]

        * **beta** (double): inverse temperature, 1/(kB*T) [ units: Ha^-1 ]

        * **Ham** (list of `nstates` lists of `nstates` doubles): diabatic state energies and couplings [ units: Ha ]


    """

    s = 208.5 * units.inv_cm2Ha # thermal energy at 300 K

    params = {}
    nosc = 25
    params["nstates"] = 3
    params["num_osc"] = nosc  # not described in the paper, but let's assume it is the same as in their 2-state model
    params["spectral_density"] = 1  # Debye
    params["model"] = 2
    params["model0"] = 2

    L, Om, T = 0.01, 0.1, 1.0

    E1 = E1 * s
    E2 = E2 * s
    V01 = V01 * s
    V02 = V02 * s
    V12 = V12 * s

    params["Omega"] = Om * s
    params["lambda"] = L * s
    params["beta"] = 1.0/(T*s)
    OM, CO, mass = GLVC.gen_bath_params(params)
    params["omega"] = [ list(OM), list(OM), list(OM) ]
    params["coupl"] = [ list(CO), list(CO), list(CO) ]
    #params["mass"] = list(mass)
    params["mass"] = [1.0 for _ in range(nosc)]

    params["Ham"] = [ [0.0,  V01, V02], 
                      [ V01, E1,  V12],
                      [ V02, V12, E2 ] ]
    params["coupling_scaling"] = [1.0, -1.0, -1.0]

    return params



def compute_model(q, params, full_id):
    model = params["model"]
    res = None

    if model==0:
        res = Esch_Levine.JCP_2020(q, params, full_id)
    elif model==1:
        res = Holstein.Holstein4(q, params, full_id)
    elif model==2:
        res = GLVC.GLVC(q, params, full_id)
    else:
        pass

    return res


def potential(q, params):
    full_id = Py2Cpp_int([0,0]) 
    
    return compute_model(q, params, full_id)


# 2-state model
model_params2_0 = {"model":0, "model0":0, "nstates":2,
                   "w0":0.015, "w1":0.005, "V":0.005, "eps":0.0, "i_crit":2, "delta":0.01 } # Esch-Levine

#model_params2_1 = {"model":0, "model0":0, "nstates":2,
#                   "w0":0.015, "w1":0.005, "V":0.005, "eps":0.0, "i_crit":2, "delta":0.01 } # Esch-Levine

# 3-state models
model_params3_0 = {"model":0, "model0":0, "nstates":3,
                   "w0":0.015, "w1":0.005, "V":0.0005, "eps":0.0, "i_crit":3, "delta":0.01 } # Esch-Levine

model_params3_1 = {"model":0, "model0":0, "nstates":3,
                   "w0":0.015, "w1":0.005, "V":0.01, "eps":0.0, "i_crit":3, "delta":0.01 } # Esch-Levine

# 5-state models
model_params5_0 = {"model":0, "model0":0, "nstates":5,
                   "w0":0.015, "w1":0.005, "V":0.005, "eps":0.0, "i_crit":4, "delta":0.01 } # Esch-Levine

model_params5_1 = {"model":0, "model0":0, "nstates":5,
                   "w0":0.015, "w1":0.005, "V":0.005, "eps":0.02, "i_crit":3, "delta":0.01 } # Esch-Levine


# Holstein models
k = 0.001
model_params6_0 = {"model":1, "model0":1, "nstates":5,
                    "E_n":[0.0, 0.0, 0.0, 0.0, 0.0 ],
                    "x_n":[0.0, 1.0, 2.0 ,3.0, 4.0],
                    "k_n":[k, 0.9*k, 0.7*k, 0.5*k, 0.25*k],
                    "V":[ [0.001, 0.000, 0.001, 0.000, 0.001 ], 
                          [0.000, 0.001, 0.000, 0.001, 0.001 ],
                          [0.001, 0.000, 0.001, 0.000, 0.001 ], 
                          [0.000, 0.001, 0.000, 0.001, 0.000 ],
                          [0.001, 0.001, 0.001, 0.000, 0.001 ]     
                       ]
                }

k = 0.001
e = -0.001
model_params6_1 = {"model":1, "model0":1, "nstates":5,
                    "E_n":[0.0, e, 2*e, 3*e, 4*e ],
                    "x_n":[0.0, 1.0, 2.0 ,3.0, 4.0],
                    "k_n":[k, k, k, k, k],
                    "V":[ [0.000, 0.000, 0.000, 0.000, 0.000 ],
                          [0.000, 0.000, 0.000, 0.000, 0.000 ],
                          [0.000, 0.000, 0.000, 0.000, 0.000 ],
                          [0.000, 0.000, 0.000, 0.000, 0.000 ],
                          [0.000, 0.000, 0.000, 0.000, 0.000 ]
                       ]
                }


# GLVS models
#                            E1   E2    V01    V02    V12 
mp0 = custom_GLVC_set(0.1, 1.0,   1.0,   1.0,   1.0);
mp1 = custom_GLVC_set(0.3, 1.0,   1.0,   1.0,   1.0);
mp2 = custom_GLVC_set(0.5, 1.0,   1.0,   1.0,   1.0);
mp3 = custom_GLVC_set(0.8, 1.0,   1.0,   1.0,   1.0);
mp4 = custom_GLVC_set(0.2, 2.0,   1.0,   1.0,   1.0);
mp5 = custom_GLVC_set(0.6, 2.0,   1.0,   1.0,   1.0);
mp6 = custom_GLVC_set(1.0, 2.0,   1.0,   1.0,   1.0);
mp7 = custom_GLVC_set(1.6, 2.0,   1.0,   1.0,   1.0);

mp0 = GLVC.get_GLVC_set4()
mp0["model"], mp0["model0"]= 2, 2

mp1 = GLVC.get_GLVC_set3()
mp1["model"], mp1["model0"] = 2, 2

all_model_params = [model_params2_0, model_params3_0, model_params3_1, model_params5_0, model_params5_1, model_params6_0, model_params6_1, 
                    mp0, mp1, mp2, mp3,    mp4, mp5, mp6, mp7  ]


