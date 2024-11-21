from mininet.net import Mininet
from mininet.link import TCLink
from mininet.topo import Topo
from sshd import sshd
from mininet.log import lg, info
import sys
import argparse

NETWORK_CONFIG = {
    "link": TCLink,
    "hosts": ['h1', 'h2', 'h3'],
    "bw": [50, 50, 50],
    "delay": ['10ms', '10ms', '10ms'],
    "ips": ['10.0.0.11', '10.0.0.12', '10.0.0.13']
}

class RoundRoleNet( Topo ):
    
    def __init__(self, *args, **params):
        Topo.__init__(self, *args, **params)
    
    def build( self ):
        self.hosts_list = []
        self.switches_list = []
        self.network_config = NETWORK_CONFIG
        
        self.switches_list.append(self.addSwitch(name='s1', stp=True, failMode='standalone'))
        
        hosts_num = len(self.network_config['hosts'])
        for i in range(hosts_num):
            h = self.addHost(name=self.network_config['hosts'][i], ip=self.network_config['ips'][i])
            self.hosts_list.append(h)
        
        # Add links, all the hosts need to link to the switch.
        for i in range(hosts_num):
            self.addLink(self.hosts_list[i], self.switches_list[0], bw=self.network_config['bw'][i], delay=self.network_config['delay'][i], use_htb=True)
        
        return        
    
def get_mininet(topo):
    return Mininet(topo, waitConnected=True, link=TCLink)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ip', type=str, help='the ip address of the host pod')
    args = parser.parse_args()
    
    lg.setLogLevel('info')
    topo = RoundRoleNet()
    net = get_mininet(topo)
    argvopts = '-D -o UseDNS=no -u0'
    
    print("argvopts: ", argvopts)
    sshd(net, opts=argvopts, routes=['10.222.0.0/24'], ip=args.ip)
