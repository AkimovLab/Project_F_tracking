import os, sys

#============= TSH ============================

for model_indx in [0, 1, 2, 3, 4, 5, 6]:
    for method_indx in [2, 4]: #[0, 1, 2, 3, 4]:
        for dt in [4.0, 20.0]:
            os.system(F"python tsh.py --model_indx={model_indx} --method_indx={method_indx} --dt={dt}")


