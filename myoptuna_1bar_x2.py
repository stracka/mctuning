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



def objective0(trial):

    mydict_temp = {        
        "top_gain":[0.1,1.5,False],        
        "top_spread":[0.001,0.3,True],        
        "bot_gain":[0.1,1.5,False],
        "bot_spread":[0.001,0.3,True],        
        "deltaE_a":[0.001,0.50,True],        
        "gain_c":[0.001,2.5,True],
        "lambda":[0.7,2,False]
    }
    
        #{
        #"mlu_scale":[0.,1.,False],
        #"c_Ej":[0.5e8,2e8,False]
        #"t_smearing":[0.1e-9,3e-9,True]
        #}

    mydict = dict()
    for key in mydict_temp.keys():
        knew = key+str(0)
        mydict.update({knew:mydict_temp[key]})
    
        #print(mydict,"\n")

    file = open("ConfigFiles/File_"+str(0)+".txt", "w")
    for key in mydict:
        val = trial.suggest_float(key,mydict[key][0],mydict[key][1],log=mydict[key][2]) #log se è true
        file.write(f"{val}\n")
    file.close()
    
    command = "root -l -b -q macro_1bar.C+\(0,1\) | grep double | awk '{print $2}' "    
    print(command)
    out = subprocess.run(command,shell=True,capture_output=True) #per  ogni command
    x = float(out.stdout.decode())
    
    return x



def objective1(trial):

    mydict_temp = {        
        "top_gain":[0.1,1.5,False],        
        "top_spread":[0.001,0.3,True],        
        "bot_gain":[0.1,1.5,False],
        "bot_spread":[0.001,0.3,True],        
        "deltaE_a":[0.001,0.50,True],        
        "gain_c":[0.001,2.5,True],
        "lambda":[0.7,2,False]
    }
    
        #{
        #"mlu_scale":[0.,1.,False],
        #"c_Ej":[0.5e8,2e8,False]
        #"t_smearing":[0.1e-9,3e-9,True]
        #}

    mydict = dict()
    for key in mydict_temp.keys():
        knew = key+str(1)
        mydict.update({knew:mydict_temp[key]})
    
        #print(mydict,"\n")

    file = open("ConfigFiles/File_"+str(1)+".txt", "w")
    for key in mydict:
        val = trial.suggest_float(key,mydict[key][0],mydict[key][1],log=mydict[key][2]) #log se è true
        file.write(f"{val}\n")
    file.close()
    
    command = "root -l -b -q macro_1bar.C+\(1,1\) | grep double | awk '{print $2}' "    
    print(command)
    out = subprocess.run(command,shell=True,capture_output=True) #per  ogni command
    x = float(out.stdout.decode())
    
    return x


def optuna_mc(n_trials=100, timeout=600):#(n_trials=500, timeout=1800): #quando fermare ottimizzazione
    """
    https://arxiv.org/pdf/1907.10902.pdf
    https://optuna.org/
    """
    
    SEED = 4005    

    logger.info("OPTUNA")

    print("hello")
    
            
    study0 = optuna.create_study(
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=SEED),
        pruner=optuna.pruners.MedianPruner(n_warmup_steps=10),
    )
    study0.optimize(objective0, n_trials=n_trials, timeout=timeout)
    #--------------------------------------------------------
    study1 = optuna.create_study(
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=SEED),
        pruner=optuna.pruners.MedianPruner(n_warmup_steps=10),
    )
    study1.optimize(objective1, n_trials=n_trials, timeout=timeout)
    
        
    #######################################################

    print("\n\nThis is the end!!!!!\n\n")

    print("BAR0_________________\n\n")

    logger.info(study0.best_trial)
    logger.info(study0.best_value)
    logger.info(study0.best_params)
    
    bestpar0 = study0.best_params 

    print("\n\nBAR1_________________\n\n")  
    logger.info(study1.best_trial)
    logger.info(study1.best_value)
    logger.info(study1.best_params)
    
    bestpar1 = study1.best_params 


    file = open("ConfigFiles/File_"+str(0)+".txt", "w")
    for key in bestpar0:
        val = bestpar0[key]
        file.write(f"{val}\n")
    file.close()


    file = open("ConfigFiles/File_"+str(1)+".txt", "w")
    for key in bestpar1:
        val = bestpar1[key]
        file.write(f"{val}\n")
    file.close()

    print("\n")
    command = "root -l -b -q macro_1bar.C+\(0,1\) | grep double | awk '{print $2}' " 
    print(command)
    subprocess.run(command,shell=True,capture_output=True)
    print("\n")

    print("\n")
    command = "root -l -b -q macro_1bar.C+\(1,1\) | grep double | awk '{print $2}' " 
    print(command)
    subprocess.run(command,shell=True,capture_output=True)
    print("\n")

    return 



#https://www.blopig.com/blog/wp-content/uploads/2019/10/GPyOpt-Tutorial1.html

#def obj_func(x):



    
#   return(out)



