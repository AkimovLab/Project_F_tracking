from liblibra_core import *
import libra_py
import libra_py.models.Esch_Levine as Esch_Levine
import libra_py.models.Holstein as Holstein


def compute_model(q, params, full_id):
    model = params["model"]
    res = None

    if model==0:
        res = Esch_Levine.JCP_2020(q, params, full_id)
    elif model==1:
        res = Holstein.Holstein4(q, params, full_id)
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


all_model_params = [model_params2_0, model_params3_0, model_params3_1, model_params5_0, model_params5_1, model_params6_0, model_params6_1   ]


