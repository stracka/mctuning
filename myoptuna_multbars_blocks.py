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

    nbars=64
    for i in range(0,nbars):
        mydict = dict()
        for key in mydict_temp.keys():
            knew = key+str(i)
            mydict.update({knew:mydict_temp[key]})
    
        #print(mydict,"\n")

        file = open("ConfigFiles/File_"+str(i)+".txt", "w")
        for key in mydict:
            val = trial.suggest_float(key,mydict[key][0],mydict[key][1],log=mydict[key][2]) #log se è true
            file.write(f"{val}\n")
        file.close()
    
    command = "root -l -b -q macro_1par_1bar.C+\(\) | grep double | awk '{print $2}' "    
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
    
            
    study = optuna.create_study(
        direction="minimize",
        sampler=optuna.samplers.TPESampler(seed=SEED),
        pruner=optuna.pruners.MedianPruner(n_warmup_steps=10),
    )
    study.optimize(objective, n_trials=n_trials, timeout=timeout)
        
    #######################################################

    print("\n\nThis is the end!!!!!\n\n")

    logger.info(study.best_trial)
    logger.info(study.best_value)
    logger.info(study.best_params)
    
    bestpar = study.best_params 

    nbars=64

    
    npars=7

    cont=0
    file = open("ConfigFiles/File_"+str(0)+".txt", "w")
    
    for key in bestpar:
        if(cont!=0 and cont%npars==0):
           file.close()
           file = open("ConfigFiles/File_"+str((int)(cont/npars))+".txt", "w")
        val = bestpar[key]
        file.write(f"{val}\n")
        cont+=1
    
    file.close()

    print("\n")
    command = "root -l -b -q macro_mlu_bar_1_2.C+\(\) | grep double | awk '{print $2}' " 
    print(command)
    subprocess.run(command,shell=True,capture_output=True)
    print("\n")

    fig = optuna.visualization.plot_intermediate_values(study) #non funzia
    fig.show()
    return study.best_trial



#https://www.blopig.com/blog/wp-content/uploads/2019/10/GPyOpt-Tutorial1.html

#def obj_func(x):



    
#   return(out)



