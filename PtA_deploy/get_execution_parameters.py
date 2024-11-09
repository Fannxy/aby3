from system_profile import *
import math

def nearest_power_of_two(n):
    # exponent = round(math.log2(n))
    exponent = math.floor(math.log2(n))
    return int(math.pow(2, exponent))


def configuration(C, M, Bp, Tb, a, n, m):
    """Compute the final coefficients.

    Args:
        C (int): # of available processord 
        M (int): overall memory limit
        Bp (int): profiler optimal block size
        Tb (int): time of evaluating ts(Bp)
        a (double): a*c + b = ti(c), corresponding coefficients.
        n (int): input sizes.
        m (int): input sizes.
    """
    # threshold = max(C * Bp, 2**20)
    # if (n*m > threshold):  # large inputs, using all the processors
    #     Bp = min((M / C), Bp, (m*n / C))
    #     return int(C), int(Bp)
    # else:  # for limited inputs, we reduce the processors.
    #     alpha = 0.2
    #     seq_eval_time = max(int(n*m / Bp), 1) * Tb * alpha
    #     c1 = int(np.sqrt(seq_eval_time / a))
    #     c1 = min(max(c1, 1), C)
    #     c1 = nearest_power_of_two(c1)
    #     return int(c1), int(min(Bp, math.ceil(n*m / c1)))    
    alpha = 0.5
    seq_eval_time = max(int(n*m / Bp), 1) * Tb * alpha
    c1 = int(math.sqrt(seq_eval_time / a))
    c1 = min(max(c1, 1), C)
    c1 = nearest_power_of_two(c1)
    return int(c1), int(min(Bp, math.ceil(n*m / c1)))

def _get_system_config():
    if(not os.path.exists(system_config_file)):
        get_deployment_configs()
    sys_cofig_dict = json.load(open(system_config_file, 'r'))
    return sys_cofig_dict['C'], sys_cofig_dict['M'], sys_cofig_dict['a'], sys_cofig_dict['b']

def _get_task_config(task):
    task_config_file = task_profiling_foler + task + "/config.json"
    
    if(not os.path.exists(task_config_file)):
        os.system(f"python {root_folder}PtA_deploy/system_profile.py --task_profile {task}")
        
    task_config_dict = json.load(open(task_config_file, 'r'))
    return task_config_dict["optimal_B"], task_config_dict["unit_time"]

def get_execution_parameters(task, n, m):
    
    C, M, a, b = _get_system_config()
    optB, Tb = _get_task_config(task)
    
    return configuration(C, M, optB, Tb, a, n, m)