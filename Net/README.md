# Setup customized network environment using mininet

## Preparation

1. Restart your pods (single one is enough) with higher privileges.

    You can modify the Deployment setting in the yaml file ('./ailab/user/userchart/templates/deployment.yaml') as:

    ```yaml
    apiVersion: apps/v1
    kind: Deployment # Add the privilege for Deployment
    metadata:
        ... # DO NOT HAVE TO CHANGE
    template:
        ... # DO NOT HAVE TO CHANGE
        containers: # Locate the `containers`
        - name: {{ $containername }}
            securityContext: # ADD THIS
            privileged: true # ADD THIS FOR HIGHER PRIVILEGE.
    ```

    Then, restart your pod by:

    ```bash
    helm upgrade release_name --values ./lab1.yaml ./userchart
    ```

    ! NOTE that all the cached code and data will be deleted, remember to store your data to PVC.

    Or directly start a new pod.


2. Start a `lo` network interface in your priviledged pod by:

    ```bash
    apt-get install -y iproute2
    ip link set lo up
    ```


3. Install mininet in your pod.

    At first install the dependencies by:

    ```bash
    cd ./mininet/
    ./util/install.sh -a
    ```

    Then, install and enable the ``openvswitch`` by:

    ```bash
    apt-get install -y openvswitch-switch
    service openvswitch-switch start
    ```

    You can test the installization by:

    ```bash
    mn --test pingall
    ```

    If the output contains something like

    ```bash
    *** Creating network
    *** Adding controller
    *** Adding hosts:
    h1 h2 
    *** Adding switches:
    s1 
    *** Adding links:
    (h1, s1) (h2, s1) 
    *** Configuring hosts
    h1 h2 
    *** Starting controller
    c0 
    *** Starting 1 switches
    s1 ...
    *** Waiting for switches to connect
    s1 
    *** Ping: testing ping reachability
    h1 -> h2 
    h2 -> h1 
    *** Results: 0% dropped (2/2 received)
    *** Stopping 1 controllers
    c0 
    *** Stopping 2 links
    ..
    *** Stopping 1 switches
    s1 
    *** Stopping 2 hosts
    h1 h2 
    *** Done
    completed in 6.554 seconds
    ```

    Then your install is success.


## Define customized ssh-based mininet environment

I provide a very simple 3-party peer-to-peer network topology through in `./mininet/examples/p2p_3pc_net.py`

The front `NETWORK_CONFIG` defines the corresponding bandwitdth (`bw`, in the unit of Mbps) and latencies. 

```python
NETWORK_CONFIG = {
    "link": TCLink,
    "hosts": ['h1', 'h2', 'h3'],
    "bw": [50, 50, 50],
    "delay": ['10ms', '10ms', '10ms'],
    "ips": ['10.0.0.11', '10.0.0.12', '10.0.0.13']
}
```

Start the network through in a `tmux` (it is necessary to keep this network environment last long because the network will shut done once the process finished.) through 

```bash
# in ./mininet/
python ./examples/p2p_3pc_net.py --ip <your pod control ip>
```

The control ip is the `10.1.0.<node-id>`.


## Run programs using mininet

The above environment will start three mininet hosts, and their ips are recorded in `./mininet/network_config.txt`.

You can ssh the mininet host through `ssh <host-ip>`. Note that the bash seems like before, you can check whether you ssh to the target host through `ifconfig`. If the network interfaces changed, it means that you are now in the new host.

(Maybe you have to wait a few seconds then can ssh to it.)

Then, we can run programs in one (main) mininet host and control the other two using `ssh`, the same as before.