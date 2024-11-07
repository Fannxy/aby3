from scipy.optimize import linprog
import itertools
import math
import numpy as np
import bisect

def calculate_expr(n, expr):
    allowed_functions = {name: getattr(math, name) for name in dir(math) if not name.startswith("__")}
    allowed_functions['n'] = n
    result = eval(expr, allowed_functions)
    return result

def calculate_mean_usage(n, mean_usage):
    keys = list(mean_usage.keys())
    keys.sort()
    pos_r = bisect.bisect_left(keys, n)
    pos_l = pos_r - 1
    if pos_l < 0:
        r = keys[pos_r]
        return n / r * mean_usage[r]
    if pos_r >= len(keys):
        l = keys[pos_l]
        return mean_usage[l]
    l = keys[pos_l]
    r = keys[pos_r]
    return (r - n) / (r - l) * mean_usage[l] + (n - l) / (r - l) * mean_usage[r]

def assign_task(strategy, bandwidth, expr_recv, expr_send, recv_mean_usage, send_mean_usage, data_size, parallelism_limit):
    num_parties = len(bandwidth)
    print(f"recv mean usage: {recv_mean_usage} | send mean usage: {send_mean_usage}")
    assert(len(expr_recv) == num_parties)
    assert(len(expr_send) == num_parties)

    coef_recv = [calculate_expr(data_size + 1, expr) - calculate_expr(data_size, expr) for expr in expr_recv]
    coef_send = [calculate_expr(data_size + 1, expr) - calculate_expr(data_size, expr) for expr in expr_send]

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

    sum_opt_x = sum(opt_x)
    print(sum_opt_x / parallelism_limit / 2)
    print(opt_x)
    indices = [i for i, x in enumerate(opt_x) if x >= sum_opt_x / parallelism_limit / 2]
    remaining_data_size = data_size
    subtask_data_size = [0] * len(perms)
    task_num = [0] * len(perms)
    for i in indices:
        subtask_data_size[i] = round(opt_x[i] / sum_opt_x * remaining_data_size)
        sum_opt_x -= opt_x[i]
        remaining_data_size -= subtask_data_size[i]
    
    predicted_bandwidth_recv = [[0] * num_parties for _ in range(len(perms))]
    predicted_bandwidth_send = [[0] * num_parties for _ in range(len(perms))]
    for party in range(num_parties):
        sum_comm_recv = 0
        sum_comm_send = 0
        for i in indices:
            sum_comm_recv += coef_recv[perms[i][party]] * subtask_data_size[i]
            sum_comm_send += coef_send[perms[i][party]] * subtask_data_size[i]
        for i in indices:
            predicted_bandwidth_recv[i][party] = coef_recv[perms[i][party]] * subtask_data_size[i] / sum_comm_recv * bandwidth[party]
            predicted_bandwidth_send[i][party] = coef_send[perms[i][party]] * subtask_data_size[i] / sum_comm_send * bandwidth[party]

    for i in indices:
        task_num[i] = parallelism_limit
        for m in range(1, parallelism_limit):
            flag = False
            for j in range(num_parties):
                predicted_mean_usage_recv = calculate_mean_usage(subtask_data_size[i] / m, recv_mean_usage[perms[i][j]])
                predicted_mean_usage_send = calculate_mean_usage(subtask_data_size[i] / m, send_mean_usage[perms[i][j]])
                # print(f"Perm {perms[i]} | Party {j} | Parallelism = {m} | Predicted mean usage (recv) = {predicted_mean_usage_recv} | Predicted mean usage (send) = {predicted_mean_usage_send} | Predicted bandwidth (recv) = {predicted_bandwidth_recv[i][j]} | Predicted bandwidth (send) = {predicted_bandwidth_send[i][j]}")
                if predicted_mean_usage_recv * m >= predicted_bandwidth_recv[i][j] or predicted_mean_usage_send * m >= predicted_bandwidth_send[i][j]:
                    flag = True
                    break
            if flag:
                task_num[i] = m
                break
    
    total_task_num = sum(task_num)
    if total_task_num > parallelism_limit:
        remaining_task_num = total_task_num
        remaining_parallelism = parallelism_limit
        for i in indices:
            tmp = round(task_num[i] / remaining_task_num * remaining_parallelism)
            remaining_task_num -= task_num[i]
            remaining_parallelism -= tmp
            task_num[i] = tmp
        total_task_num = sum(task_num)

    # for i in indices:
    #     print(f"Task {perms[i]}: ")
    #     print(f"    Data size: {subtask_data_size[i]}")
    #     print(f"    Task num: {task_num[i]}")
    #     print(f"    Predicted bandwidth (recv): {predicted_bandwidth_recv[i]}")
    #     print(f"    Predicted bandwidth (send): {predicted_bandwidth_send[i]}")

    result = []
    if strategy == 'baseline':
        result = [[perms[0], data_size // total_task_num] for _ in range(total_task_num)]
    elif strategy == 'roundrole':
        for i in indices:
            for _ in range(task_num[i]):
                result.append([perms[i], data_size // total_task_num])
    else:
        raise ValueError(f"Invalid strategy: {strategy}")

    return result
