#!/usr/bin/env python
# coding: utf-8

"""
Optimization with various strategies: 
For illustration, currently limited to tree ensemble classifiers

- Optuna: tree-structured Parzen estimator (TPE). Not limited to sklearn

"""

    
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import json
import subprocess

import logging
FORMAT = '%(asctime)-15s- %(levelname)s - %(name)s -%(message)s'
logging.basicConfig(format=FORMAT, level=logging.INFO)
logger = logging.getLogger(__name__)


# optuna
import optuna

# gpyopt
#import GPyOpt
#from GPyOpt.methods import BayesianOptimization



def objective(trial):

    mydict = {        
        "mlu_scale":[0.,1.,False],
        "top_gain":[0.1,1.5,False],        
        "top_spread":[0.001,0.3,True],        
        "bot_gain":[0.1,1.5,False],
        "bot_spread":[0.001,0.3,True],        
        "deltaE_a":[0.001,0.50,True],        
        "gain_c":[0.001,2.5,True],
        "lambda":[0.7,2,False],
    }

    params=''
    for key in mydict:
        val = trial.suggest_float(key,mydict[key][0],mydict[key][1],log=mydict[key][2])
        params += str(val) + ','
    params=params[:-1]

    command = "root -l -b -q macro.C+\("+ params +"\) | grep double | awk '{print $2}' " 

    print(command)
    out = subprocess.run(command,shell=True,capture_output=True)
    x = float(out.stdout.decode())
    
    return x



def optuna_mc(n_trials=100, timeout=600): 
    """
    https://arxiv.org/pdf/1907.10902.pdf
    https://optuna.org/
    """
    
    SEED = 4005    

    logger.info("OPTUNA")
    
            
    study = optuna.create_study(
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=SEED),
        pruner=optuna.pruners.MedianPruner(n_warmup_steps=10),
    )
    study.optimize(objective, n_trials=n_trials, timeout=timeout)
        
    logger.info(study.best_trial)
    logger.info(study.best_value)
    logger.info(study.best_params)
    
    return study.best_trial



#https://www.blopig.com/blog/wp-content/uploads/2019/10/GPyOpt-Tutorial1.html

#def obj_func(x):



    
#   return(out)



