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

def assign_aggr(assignment):
    if len(assignment) == 1:
        return [assignment]
    left_assignment = assign_aggr(assignment[:len(assignment) // 2])
    right_assignment = assign_aggr(assignment[len(assignment) // 2:])
    new_assignment = []
    for i in range(max(len(left_assignment), len(right_assignment))):
        new_layer = []
        if i < len(left_assignment):
            new_layer += left_assignment[i]
        if i < len(right_assignment):
            new_layer += right_assignment[i]
        new_assignment.append(new_layer)
    new_assignment.append([[assignment[0][0], sum([item[1] for item in assignment])]])
    return new_assignment

def assign_task(strategy, bandwidth, expr_recv, expr_send, recv_mean_usage, send_mean_usage, data_size, parallelism_limit, agg_time):
    num_parties = len(bandwidth)
    # print(f"recv mean usage: {recv_mean_usage} | send mean usage: {send_mean_usage}")
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
    indices = [i for i, x in enumerate(opt_x) if x >= sum_opt_x / parallelism_limit / 2]
    remaining_data_size = data_size
    subtask_data_size = [0] * len(perms)
    task_num = [0] * len(perms)
    for i in indices:
        subtask_data_size[i] = round(opt_x[i] / sum_opt_x * remaining_data_size)
        sum_opt_x -= opt_x[i]
        remaining_data_size -= subtask_data_size[i]
    
    min_time = 1e9
    # for m in range(1, parallelism_limit + 1):
    for m in range(parallelism_limit, parallelism_limit + 1):
        m_task_num = [0] * len(perms)
        remaining_data_size = data_size
        remaining_parallelism = m
        for i in indices:
            m_task_num[i] = round(subtask_data_size[i] / remaining_data_size * remaining_parallelism)
            remaining_data_size -= subtask_data_size[i]
            remaining_parallelism -= m_task_num[i]
        
        sum_comm_recv = [0] * num_parties
        sum_comm_send = [0] * num_parties
        sum_mean_usage_recv = [0] * num_parties
        sum_mean_usage_send = [0] * num_parties
        for i in indices:
            if m_task_num[i] == 0:
                continue
            for j in range(num_parties):
                sum_comm_recv[j] += coef_recv[perms[i][j]] * data_size / m * m_task_num[i]
                sum_comm_send[j] += coef_send[perms[i][j]] * data_size / m * m_task_num[i]
                sum_mean_usage_recv[j] += calculate_mean_usage(data_size / m, recv_mean_usage[perms[i][j]]) / (2 * bandwidth[j]) * m_task_num[i]
                sum_mean_usage_send[j] += calculate_mean_usage(data_size / m, send_mean_usage[perms[i][j]]) / (2 * bandwidth[j]) * m_task_num[i]
        
        max_usage = 1
        for j in range(num_parties):
            max_usage = max(max_usage, sum_mean_usage_recv[j], sum_mean_usage_send[j])
        for j in range(num_parties):
            sum_mean_usage_recv[j] /= max_usage
            sum_mean_usage_send[j] /= max_usage
        
        time = 0
        for j in range(num_parties):
            if sum_comm_recv[j] > 0:
                time = max(time, sum_comm_recv[j] / (bandwidth[j] * sum_mean_usage_recv[j]))
            if sum_comm_send[j] > 0:
                time = max(time, sum_comm_send[j] / (bandwidth[j] * sum_mean_usage_send[j]))
        
        print(f"m = {m}")
        print(f"    computation time = {time}")
        for j in range(num_parties):
            print(f"    Party {j}: recv = {sum_comm_recv[j]} | send = {sum_comm_send[j]} | recv usage = {sum_mean_usage_recv[j]} | send usage = {sum_mean_usage_send[j]}")

        agg_task_size = m
        while agg_task_size > 1:
            agg_task_size = math.ceil(agg_task_size / 2)
            time += calculate_mean_usage(agg_task_size * data_size / m, agg_time)
        
        print(f"    total time = {time}")
        
        if time < min_time:
            min_time = time
            task_num = m_task_num

    total_task_num = sum(task_num)

    # for i in indices:
    #     print(f"Task {perms[i]}: ")
    #     print(f"    Data size: {subtask_data_size[i]}")
    #     print(f"    Task num: {task_num[i]}")

    comp_assignment = []
    aggr_assignment = []
    if strategy == 'baseline':
        comp_assignment = [[perms[0], data_size // total_task_num] for _ in range(total_task_num)]
    elif strategy == 'roundrole':
        for i in indices:
            for _ in range(task_num[i]):
                comp_assignment.append([perms[i], data_size // total_task_num])
    else:
        raise ValueError(f"Invalid strategy: {strategy}")
    
    # aggr_layer = comp_assignment
    # while len(aggr_layer) > 1:
    #     new_layer = []
    #     for i in range(0, len(aggr_layer) // 2):
    #         new_layer.append([aggr_layer[i][0], aggr_layer[i][1] + aggr_layer[len(aggr_layer) // 2 * 2 - i - 1][1]])
    #     aggr_assignment.append([item for item in new_layer])
    #     if len(aggr_layer) % 2 == 1:
    #         new_layer.append([aggr_layer[-1][0], aggr_layer[-1][1]])
    #     new_layer.sort(key=lambda x: x[1])
    #     aggr_layer = new_layer

    aggr_assignment = assign_aggr(comp_assignment)
    
    return comp_assignment, aggr_assignment
