############################################################################
# Copyright (c) 2014-2018 Beijing Institutes of Life Science
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import numpy as np
import math
import copy

def logistic_func(theta, x):
    return float(1) / (1 + math.e**(-x.dot(theta)))

def log_gradient(theta, x, y):
    first_calc = logistic_func(theta, x) - np.squeeze(y)
    final_calc = first_calc.T.dot(x)
    return final_calc

def cost_func(theta, x, y):
    log_func_v = logistic_func(theta,x)
    y = np.squeeze(y)
    step1 = y * np.log(log_func_v)
    step2 = (1-y) * np.log(1 - log_func_v)
    final = -step1 - step2
    return np.mean(final)

def grad_desc(theta_values, X, y, lr=.001, converge_change=.001):
    #normalize
    X = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
    #setup cost iter
    cost_iter = []
    cost = cost_func(theta_values, X, y)
    cost_iter.append([0, cost])
    change_cost = 1
    i = 1
    while(change_cost > converge_change):
        old_cost = cost
        theta_values = theta_values - (lr * log_gradient(theta_values, X, y))
        cost = cost_func(theta_values, X, y)
        cost_iter.append([i, cost])
        change_cost = old_cost - cost
        i+=1
    return theta_values, np.array(cost_iter)

def pred_values(theta, X, hard=True):
    #normalize
    X = (X - np.mean(X, axis=0)) / np.std(X, axis=0)
    pred_prob = logistic_func(theta, X)
    pred_value = np.where(pred_prob >= .45, 1, 0)
    if hard:
        return pred_value
    return pred_prob

def outliers(obs):
    tmp = copy.copy(obs)
    for i in range(5):
        std = np.array(tmp).std()
        mean = np.array(tmp).mean()
        tmp = [a for a in tmp if mean-3*std <= a <=mean+3*std]
    return round(np.array(tmp).mean(),4), round(np.array(tmp).std(),4)

flag = {}
with open('totalBubLenRate') as af:
    for l in af:
        ae = l.rstrip().split()
        if ae[0] not in flag :
            flag[ae[0]] = []
        flag[ae[0]].append(float(ae[1]))

truth = {}
with open('flags') as af:
    for line in af:
        ae = line.rstrip().split()
        truth[ae[0]] = int(ae[-1])


tx, ty = [], []
for ele in flag:
    m, t =  outliers(flag[ele])
    tx.append([m, round(math.log(m/len(flag[ele])), 4)])
    ty.append(truth[ele])

train_x = np.array(tx)
train_y = np.array(ty)

betas = np.zeros(train_x.shape[1])
# training to get betas
fitted_values, cost_iter = grad_desc(betas, train_x, train_y)
# test
predicted_y = pred_values(fitted_values, train_x)
# accuracy
np.sum(train_y == predicted_y)

