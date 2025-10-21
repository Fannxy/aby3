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

    Install mininet by:

    ```bash
    git clone git@github.com:Fannxy/mininet.git
    cp ./install_wo_sudo.sh ./mininet/util/
    ```
    
    At first install the dependencies by:

    ```bash
    cd ./mininet/
    ./util/install.sh -a
    ```
    (You can use `tmux` since this will take a while.)
   
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


    Install the requirement from opensource version: https://docs.openvswitch.org/en/latest/intro/install/general/

If failed! see `## Problems with BTF Errors`, which records some of my previous attempts,  and ask for GPT with help :) 233333. This problem is highly correlated with your physical nodes so...

## Define customized ssh-based mininet environment

I provide a very simple 3-party peer-to-peer network topology through in `p2p_3pc_net.py`, copy it into the `./mininet/examples/`

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


## Problems with BTF Errors

Currently, the 5.15.0-127-generic version linux kernel have some problems with the nf_conncount module. 

```
> dmseg | tail
failed to validate module [nf_conncount] BTF: -22
```

In this case, we can manually start the ovs vswitch to bypass the BTF problems.

1. Install the ovs through the source code, following: https://docs.openvswitch.org/en/latest/intro/install/general/

!!! Note that the version currently being tests is v3.2.0! Then follow the commands in the previous doc to build from source till `make install`.

2. Try to start `ovs` through the Starting recommendations in the previous doc. If you succeed, congratulations!

You can use `mn --test pingall` to test whether you succeed or not!


4. If failed, we can also use the following commands to maunally start the vswitch (succeeded once, but I do not know why...) In this case, ask GPT for detailed help!

    创建数据库目录
    mkdir -p /usr/local/etc/openvswitch
    mkdir -p /usr/local/var/run/openvswitch

    初始化 OVS 数据库
    ovsdb-tool create /usr/local/etc/openvswitch/conf.db vswitchd/vswitch.ovsschema

    启动 ovsdb-server
    ovsdb-server /usr/local/etc/openvswitch/conf.db --remote=punix:/usr/local/var/run/openvswitch/db.sock --remote=db:Open_vSwitch,Open_vSwitch,manager_options --pidfile=/usr/local/var/run/openvswitch/ovsdb-server.pid --detach --log-file=/var/log/openvswitch/ovsdb-server.log

    初始化 OVS 控制数据库（只需执行一次）
    ovs-vsctl --db=unix:/usr/local/var/run/openvswitch/db.sock --no-wait init

    启动 ovs-vswitchd（使用 netdev 用户态 datapath）
    ovs-vswitchd --pidfile --detach --log-file=/var/log/openvswitch/ovs-vswitchd.log --unixctl=/usr/local/var/run/openvswitch/vswitchd.sock

5. In this case, using the following code to define the virtual switch.

    ```
    s1 = net.addSwitch('s1', datapath='user', cls=OVSSwitch)
    ```
