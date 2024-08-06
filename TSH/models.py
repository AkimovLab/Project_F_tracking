from liblibra_core import *
import libra_py
import libra_py.models.Esch_Levine as Esch_Levine


def compute_model(q, params, full_id):
    model = params["model"]
    res = None

    if model==0:
        res = Esch_Levine.JCP_2020(q, params, full_id)
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


all_model_params = [model_params2_0, model_params3_0, model_params3_1, model_params5_0, model_params5_1  ]


