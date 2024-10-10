from scipy.optimize import linprog
import itertools
import math
import numpy as np

def calculate_expression(n, expression):
    allowed_functions = {name: getattr(math, name) for name in dir(math) if not name.startswith("__")}
    allowed_functions['n'] = n
    result = eval(expression, allowed_functions)
    return result

def assign_task(bandwidth, expr_recv, expr_send, data_size):
    num_parties = len(bandwidth)
    assert(len(expr_recv) == num_parties)
    assert(len(expr_send) == num_parties)

    coef_recv = [calculate_expression(data_size + 1, expr) - calculate_expression(data_size, expr) for expr in expr_recv]
    coef_send = [calculate_expression(data_size + 1, expr) - calculate_expression(data_size, expr) for expr in expr_send]

    perms = list(itertools.permutations([i for i in range(num_parties)])) # role of i-th party is perm[i]

    c = np.zeros(len(perms) + 1)
    c[-1] = 1

    A_ub = []
    for i in range(num_parties):
        u_recv = []
        for perm in perms:
            u_recv.append(coef_recv[perm[i]])
        u_recv.append(-bandwidth[i])
        A_ub.append(u_recv)

        u_send = []
        for perm in perms:
            u_send.append(coef_send[perm[i]])
        u_send.append(-bandwidth[i])
        A_ub.append(u_send)
    A_ub = np.array(A_ub)
    b_ub = np.zeros(num_parties * 2)

    A_eq = np.zeros((1, len(perms) + 1))
    A_eq[0, :-1] = 1
    b_eq = np.array([data_size])
    
    bounds = [(0, None)] * (len(perms)+1)
    
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
    opt_t = res.x[-1]
    opt_x = res.x[:-1]

    threshold = max(opt_x) / num_parties
    indices = [i for i, x in enumerate(opt_x) if x >= threshold]
    sum_opt_x = sum(opt_x)
    remaining_data_size = data_size
    result = []
    for i in indices:
        subtask_data_size = round(opt_x[i] / sum_opt_x * remaining_data_size)
        sum_opt_x -= opt_x[i]
        remaining_data_size -= subtask_data_size
        result.append([perms[i], subtask_data_size])
    return result
