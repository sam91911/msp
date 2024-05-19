import numpy as np
import copy
from tqdm import tqdm

class sp:
    def __init__(self):
        self.r = {}
        self.c = {}
        self.l = {}
        self.vsrc = {}
        self.isrc = {}
        self.esrc = {}
        self.gsrc = {}
        self.etype = {}
        self.f = []
        self.nodes = ['0']
        self.nodes_id = {'0': 0}
        self.cnt = 0
    def R(self, src, dst, value, name=None):
        if name is None:
            name = str(self.cnt)
            self.cnt += 1
        if src not in self.nodes:
            self.nodes_id[src] = len(self.nodes)
            self.nodes.append(src)
        if dst not in self.nodes:
            self.nodes_id[dst] = len(self.nodes)
            self.nodes.append(dst)
        self.r[name] = [self.nodes_id[src], self.nodes_id[dst], value]
        self.etype[name] = "R"
        return
    def C(self, src, dst, value, ic=0, name=None):
        if name is None:
            name = str(self.cnt)
            self.cnt += 1
        if src not in self.nodes:
            self.nodes_id[src] = len(self.nodes)
            self.nodes.append(src)
        if dst not in self.nodes:
            self.nodes_id[dst] = len(self.nodes)
            self.nodes.append(dst)
        self.c[name] = [self.nodes_id[src], self.nodes_id[dst], value, ic]
        self.etype[name] = "C"
        return
    def L(self, src, dst, value, ic=0, name=None):
        if name is None:
            name = str(self.cnt)
            self.cnt += 1
        if src not in self.nodes:
            self.nodes_id[src] = len(self.nodes)
            self.nodes.append(src)
        if dst not in self.nodes:
            self.nodes_id[dst] = len(self.nodes)
            self.nodes.append(dst)
        self.l[name] = [self.nodes_id[src], self.nodes_id[dst], value, ic]
        self.etype[name] = "L"
        return
    def Vsrc(self, src, dst, value, name=None):
        if name is None:
            name = str(self.cnt)
            self.cnt += 1
        if src not in self.nodes:
            self.nodes_id[src] = len(self.nodes)
            self.nodes.append(src)
        if dst not in self.nodes:
            self.nodes_id[dst] = len(self.nodes)
            self.nodes.append(dst)
        self.vsrc[name] = [self.nodes_id[src], self.nodes_id[dst], value]
        self.etype[name] = "V"
        return
    def Isrc(self, src, dst, value, name=None):
        if name is None:
            name = str(self.cnt)
            self.cnt += 1
        if src not in self.nodes:
            self.nodes_id[src] = len(self.nodes)
            self.nodes.append(src)
        if dst not in self.nodes:
            self.nodes_id[dst] = len(self.nodes)
            self.nodes.append(dst)
        self.isrc[name] = [self.nodes_id[src], self.nodes_id[dst], value]
        self.etype[name] = "I"
        return
    def Esrc(self, src, dst, value, np, nn, name=None):
        if name is None:
            name = str(self.cnt)
            self.cnt += 1
        if src not in self.nodes:
            self.nodes_id[src] = len(self.nodes)
            self.nodes.append(src)
        if dst not in self.nodes:
            self.nodes_id[dst] = len(self.nodes)
            self.nodes.append(dst)
        if np not in self.nodes:
            self.nodes_id[np] = len(self.nodes)
            self.nodes.append(np)
        if nn not in self.nodes:
            self.nodes_id[nn] = len(self.nodes)
            self.nodes.append(nn)
        self.esrc[name] = [self.nodes_id[src], self.nodes_id[dst], value, self.nodes_id[np], self.nodes_id[nn]]
        self.etype[name] = "E"
        return
    def Gsrc(self, src, dst, value, np, nn, name=None):
        if name is None:
            name = str(self.cnt)
            self.cnt += 1
        if src not in self.nodes:
            self.nodes_id[src] = len(self.nodes)
            self.nodes.append(src)
        if dst not in self.nodes:
            self.nodes_id[dst] = len(self.nodes)
            self.nodes.append(dst)
        if np not in self.nodes:
            self.nodes_id[np] = len(self.nodes)
            self.nodes.append(np)
        if nn not in self.nodes:
            self.nodes_id[nn] = len(self.nodes)
            self.nodes.append(nn)
        self.gsrc[name] = [self.nodes_id[src], self.nodes_id[dst], value, self.nodes_id[np], self.nodes_id[nn]]
        self.etype[name] = "G"
        return
    def add_f(self, f):
        self.f.append(f)
        return
    def root_node(self, node):
        if self.cnodes[node] == node:
            return node
        else:
            self.cnodes[node] = self.root_node(self.cnodes[node])
            return self.cnodes[node]
    def root_merge(self, node1, node2):
        node1 = self.root_node(node1)
        node2 = self.root_node(node2)
        if node1 == 0:
            u = 0
            d = node2
        elif node2 == 0:
            u = 0
            d = node1
        elif self.rank[node1] >= self.rank[node2]:
            u = node1
            d = node2
        else:
            u = node2
            d = node1
        if(u == d):
            raise "voltage source loop"
            return
        if self.rank[u] == self.rank[d]:
            self.rank[u] += 1
        self.cnodes[d] = u
        return
    def setup(self):
        print("start setup")
        self.cnodes = list(range(len(self.nodes)))
        self.rank = [0]*len(self.nodes)
        nodemap = {}
        self.nmap = []
        cnt = 0
        for x in self.c.values():
            self.root_merge(x[0], x[1])
        for x in self.vsrc.values():
            self.root_merge(x[0], x[1])
        for x in self.esrc.values():
            self.root_merge(x[0], x[1])
        for x in range(len(self.nodes)):
            self.root_node(x)
        for x in self.cnodes:
            if x not in nodemap:
                nodemap[x] = cnt
                cnt += 1
        self.nlen = cnt
        for x in self.cnodes:
            self.nmap.append(nodemap[self.root_node(x)])
        self.vdep = {}
        for i in range(len(self.nodes)):
            self.vdep[i] = set()
        depend = list([i in self.cnodes for i in range(len(self.nodes))])
        while (not all(depend)):
            for x in self.c.values():
                src = x[0]
                dst = x[1]
                if (depend[src])and(not depend[dst]):
                    depend[dst] = True
                    self.vdep[src].add(dst)
                if (depend[dst])and(not depend[src]):
                    depend[src] = True
                    self.vdep[dst].add(src)
            for x in self.vsrc.values():
                src = x[0]
                dst = x[1]
                if (depend[src])and(not depend[dst]):
                    depend[dst] = True
                    self.vdep[src].add(dst)
                if (depend[dst])and(not depend[src]):
                    depend[src] = True
                    self.vdep[dst].add(src)
            for x in self.esrc.values():
                src = x[0]
                dst = x[1]
                if (depend[src])and(not depend[dst]):
                    depend[dst] = True
                    self.vdep[src].add(dst)
                if (depend[dst])and(not depend[src]):
                    depend[src] = True
                    self.vdep[dst].add(src)
        for x in self.vdep:
            depend = copy.copy(self.vdep[x])
            while(len(depend) != 0):
                temp = depend.pop()
                depend |= self.vdep[temp]
                self.vdep[x] |= self.vdep[temp]
        self.fG = np.zeros((len(self.nodes), len(self.nodes)))
        for x in self.r.values():
            self.fG[x[0], x[0]] += 1/x[2]
            self.fG[x[0], x[1]] -= 1/x[2]
            self.fG[x[1], x[1]] += 1/x[2]
            self.fG[x[1], x[0]] -= 1/x[2]
        for x in self.gsrc.values():
            self.fG[x[0], x[3]] -= x[2]
            self.fG[x[1], x[3]] += x[2]
            self.fG[x[0], x[4]] += x[2]
            self.fG[x[1], x[4]] -= x[2]
        self.sG = np.zeros((self.nlen, len(self.nodes)))
        for x in self.r.values():
            self.sG[self.nmap[x[0]], x[0]] += 1/x[2]
            self.sG[self.nmap[x[0]], x[1]] -= 1/x[2]
            self.sG[self.nmap[x[1]], x[1]] += 1/x[2]
            self.sG[self.nmap[x[1]], x[0]] -= 1/x[2]
        for x in self.gsrc.values():
            self.sG[self.nmap[x[0]], x[3]] -= x[2]
            self.sG[self.nmap[x[1]], x[3]] += x[2]
            self.sG[self.nmap[x[0]], x[4]] += x[2]
            self.sG[self.nmap[x[1]], x[4]] -= x[2]
        self.G = np.zeros((self.nlen, self.nlen))
        for x in self.r.values():
            self.G[self.nmap[x[0]], self.nmap[x[0]]] += 1/x[2]
            self.G[self.nmap[x[0]], self.nmap[x[1]]] -= 1/x[2]
            self.G[self.nmap[x[1]], self.nmap[x[1]]] += 1/x[2]
            self.G[self.nmap[x[1]], self.nmap[x[0]]] -= 1/x[2]
        for x in self.gsrc.values():
            self.G[self.nmap[x[0]], self.nmap[x[3]]] -= x[2]
            self.G[self.nmap[x[1]], self.nmap[x[3]]] += x[2]
            self.G[self.nmap[x[0]], self.nmap[x[4]]] += x[2]
            self.G[self.nmap[x[1]], self.nmap[x[4]]] -= x[2]
        self.videp = []
        depend = [0]*len(self.nodes)
        for x in self.c.values():
            depend[x[0]] += 1
            depend[x[1]] += 1
        for x in self.vsrc.values():
            depend[x[0]] += 1
            depend[x[1]] += 1
        for x in self.esrc.values():
            depend[x[0]] += 1
            depend[x[1]] += 1
        while (any(depend)):
            for x in self.c:
                src = self.c[x][0]
                dst = self.c[x][1]
                if x not in self.videp:
                    if depend[src] == 1:
                        depend[src] -= 1
                        depend[dst] -= 1
                        self.videp.append(("C", x, True))
                    elif depend[dst] == 1:
                        depend[src] -= 1
                        depend[dst] -= 1
                        self.videp.append(("C", x, False))
            for x in self.vsrc:
                src = self.vsrc[x][0]
                dst = self.vsrc[x][1]
                if x not in self.videp:
                    if depend[src] == 1:
                        depend[src] -= 1
                        depend[dst] -= 1
                        self.videp.append(("V", x, True))
                    elif depend[dst] == 1:
                        depend[src] -= 1
                        depend[dst] -= 1
                        self.videp.append(("V", x, False))
            for x in self.esrc:
                src = self.esrc[x][0]
                dst = self.esrc[x][1]
                if x not in self.videp:
                    if depend[src] == 1:
                        depend[src] -= 1
                        depend[dst] -= 1
                        self.videp.append(("E", x, True))
                    elif depend[dst] == 1:
                        depend[src] -= 1
                        depend[dst] -= 1
                        self.videp.append(("E", x, False))
        depend = {}
        self.edep = []
        for x in self.esrc:
            depend[x] = []
            src = self.esrc[x][0]
            dst = self.esrc[x][1]
            for y in self.esrc:
                dep = self.esrc[y][0]
                if (src in self.vdep[dep])^(dst in self.vdep[dep]):
                    depend[x].append[y]
        while(len(self.edep) != len(self.esrc)):
            act = False
            for x in self.esrc:
                if(len(depend[x]) == 0)and(x not in self.edep):
                    self.edep.append(x)
                    act = True
                    for y in self.esrc:
                        if(x in depend[y]):
                            depend[y].remove(x)
            if not act:
                raise "dependent voltage source loop"
        self.Ae = np.identity(len(self.nodes))
        for x in self.edep:
            it = self.esrc[x]
            if it[0] in self.vdep[it[1]]:
                self.Ae[it[0], ::] += it[2]*self.Ae[it[3], ::]
                self.Ae[it[0], ::] -= it[2]*self.Ae[it[4], ::]
                for y in self.vdep[it[0]]:
                    self.Ae[y, ::] += it[2]*self.Ae[it[3], ::]
                    self.Ae[y, ::] -= it[2]*self.Ae[it[4], ::]
            else:
                self.Ae[it[1], ::] -= it[2]*self.Ae[it[3], ::]
                self.Ae[it[1], ::] += it[2]*self.Ae[it[4], ::]
                for y in self.vdep[it[1]]:
                    self.Ae[y, ::] -= it[2]*self.Ae[it[3], ::]
                    self.Ae[y, ::] += it[2]*self.Ae[it[4], ::]
        self.Ve = np.zeros((len(self.nodes), self.nlen))
        for x in self.esrc.values():
            if x[1] in self.vdep[x[0]]:
                self.Ve[x[1], self.nmap[x[3]]] -= x[2]
                self.Ve[x[1], self.nmap[x[4]]] += x[2]
                for y in self.vdep[x[1]]:
                    self.Ve[y, self.nmap[x[3]]] -= x[2]
                    self.Ve[y, self.nmap[x[4]]] += x[2]
            else:
                self.Ve[x[0], self.nmap[x[3]]] += x[2]
                self.Ve[x[0], self.nmap[x[4]]] -= x[2]
                for y in self.vdep[x[0]]:
                    self.Ve[y, self.nmap[x[3]]] += x[2]
                    self.Ve[y, self.nmap[x[4]]] -= x[2]
        self.Iid = np.zeros((self.nlen,))
        for x in self.l.values():
            self.Iid[self.nmap[x[0]]] -= x[3]
            self.Iid[self.nmap[x[1]]] += x[3]
        for x in self.isrc.values():
            self.Iid[self.nmap[x[0]]] -= x[2]
            self.Iid[self.nmap[x[1]]] += x[2]
        self.Vid = np.zeros((len(self.nodes),))
        for x in self.c.values():
            if x[0] in self.vdep[x[1]]:
                self.Vid[x[0]] += x[3]
                for y in self.vdep[x[0]]:
                    self.Vid[y] += x[3]
            else:
                self.Vid[x[1]] -= x[3]
                for y in self.vdep[x[1]]:
                    self.Vid[y] -= x[3]
        for x in self.vsrc.values():
            if x[0] in self.vdep[x[1]]:
                self.Vid[x[0]] += x[2]
                for y in self.vdep[x[0]]:
                    self.Vid[y] += x[2]
            else:
                self.Vid[x[1]] -= x[2]
                for y in self.vdep[x[1]]:
                    self.Vid[y] -= x[2]
        print("end setup")
        return
    def get_V(self, update = True):
        G = (self.G+(self.sG@self.Ve))[1::, 1::]
        I = (self.Iid-(self.sG@(self.Ae@self.Vid)))[1::]
        V = np.linalg.lstsq(G, I, rcond = None)[0]
        vnet = np.copy(self.Vid)
        for i in range(len(vnet)):
            if self.nmap[i] == 0:
                continue
            vnet[i] += V[self.nmap[i]-1]
        vnet = self.Ae@vnet
        if update:
            self.vnet = vnet
        return vnet
    def get_I(self, update = True, V = None):
        I = np.zeros((len(self.nodes),))
        for x in self.l.values():
            I[x[0]] -= x[3]
            I[x[1]] += x[3]
        for x in self.isrc.values():
            I[x[0]] -= x[2]
            I[x[1]] += x[2]
        if V is None:
            V = self.get_V()
        I -= self.fG@V
        V_I = {}
        for x in self.videp:
            if x[0] == "C":
                src = self.c[x[1]][0]
                dst = self.c[x[1]][1]
                if x[2]:
                    ei = I[src]
                else:
                    ei = -I[dst]
                V_I[x[1]] = ei
                I[src] -= ei
                I[dst] += ei
            if x[0] == "V":
                src = self.vsrc[x[1]][0]
                dst = self.vsrc[x[1]][1]
                if x[2]:
                    ei = I[src]
                else:
                    ei = -I[dst]
                V_I[x[1]] = ei
                I[src] -= ei
                I[dst] += ei
            if x[0] == "E":
                src = self.esrc[x[1]][0]
                dst = self.esrc[x[1]][1]
                if x[2]:
                    ei = I[src]
                else:
                    ei = -I[dst]
                V_I[x[1]] = ei
                I[src] -= ei
                I[dst] += ei
        if update:
            self.V_I = V_I
        return V_I
    def run(self):
        self.get_I()
        return
    def update_d(self, I, V):
        G = (self.G+(self.sG@self.Ve))[1::, 1::]
        I = (I-(self.sG@(self.Ae@V)))[1::]
        Vans = np.linalg.lstsq(G, I, rcond = None)[0]
        vnet = np.copy(self.Vid)
        for i in range(len(vnet)):
            if self.nmap[i] == 0:
                continue
            vnet[i] += Vans[self.nmap[i]]
        vnet = self.Ae@vnet
        I = np.zeros((len(self.nodes),))
        for x in self.l.values():
            I[x[0]] -= x[3]
            I[x[1]] += x[3]
        for x in self.isrc.values():
            I[x[0]] -= x[2]
            I[x[1]] += x[2]
        V = vnet
        I -= self.fG@V
        V_I = {}
        for x in self.videp:
            if x[0] == "C":
                src = self.c[x[1]][0]
                dst = self.c[x[1]][1]
                if x[2]:
                    ei = I[src]
                else:
                    ei = -I[dst]
                V_I[x[1]] = ei
                I[src] -= ei
                I[dst] += ei
            if x[0] == "V":
                src = self.vsrc[x[1]][0]
                dst = self.vsrc[x[1]][1]
                if x[2]:
                    ei = I[src]
                else:
                    ei = -I[dst]
                V_I[x[1]] = ei
                I[src] -= ei
                I[dst] += ei
            if x[0] == "E":
                src = self.esrc[x[1]][0]
                dst = self.esrc[x[1]][1]
                if x[2]:
                    ei = I[src]
                else:
                    ei = -I[dst]
                V_I[x[1]] = ei
                I[src] -= ei
                I[dst] += ei
        C_dv = {}
        L_di = {}
        for x in self.c:
            C_dv[x] = V_I[x]/self.c[x][2]
        for x in self.l:
            dv = V[self.l[x][0]] - V[self.l[x][1]]
            L_di[x] = dv/self.l[x][2]
        return C_dv, L_di
    def update_s(self):
        V = self.vnet
        V_I= self.V_I
        C_dv = {}
        L_di = {}
        for x in self.c:
            C_dv[x] = V_I[x]/self.c[x][2]
        for x in self.l:
            dv = V[self.l[x][0]] - V[self.l[x][1]]
            L_di[x] = dv/self.l[x][2]
        return C_dv, L_di
    def update(self, dt, t, C_dv, L_di):
        for x in self.c:
            src = self.c[x][0]
            dst = self.c[x][1]
            dV = C_dv[x]*dt
            if src in self.vdep[dst]:
                self.Vid[src] += dV
                for y in self.vdep[src]:
                    self.Vid[y] += dV
            else:
                self.Vid[dst] -= dV
                for y in self.vdep[dst]:
                    self.Vid[y] -= dV
            self.c[x][3] += C_dv[x]*dt
        for x in self.l:
            src = self.l[x][0]
            dst = self.l[x][1]
            dI = L_di[x]*dt
            self.Iid[self.nmap[src]] -= dI
            self.Iid[self.nmap[dst]] += dI
            self.l[x][3] += L_di[x]*dt
        for f in self.f:
            f(self, dt, t)
        return
    def deepcopy(self, dst = None):
        if dst is None:
            dst = sp()
        dst.r = copy.deepcopy(self.r)
        dst.c = copy.deepcopy(self.c)
        dst.l = copy.deepcopy(self.l)
        dst.isrc = copy.deepcopy(self.isrc)
        dst.vsrc = copy.deepcopy(self.vsrc)
        dst.f = copy.deepcopy(self.f)
        dst.nodes = copy.deepcopy(self.nodes)
        dst.nodes_id = copy.deepcopy(self.nodes_id)
        dst.vdep = copy.deepcopy(self.vdep)
        dst.videp = copy.deepcopy(self.videp)
        dst.edep = copy.deepcopy(self.edep)
        dst.cnodes = copy.deepcopy(self.cnodes)
        dst.rank = copy.deepcopy(self.rank)
        dst.nmap = copy.deepcopy(self.nmap)
        dst.nlen = copy.deepcopy(self.nlen)
        dst.cnt = copy.deepcopy(self.cnt)
        dst.esrc = copy.deepcopy(self.esrc)
        dst.gsrc = copy.deepcopy(self.gsrc)
        dst.fG = np.copy(self.fG)
        dst.sG = np.copy(self.sG)
        dst.G = np.copy(self.G)
        dst.Ve = np.copy(self.Ve)
        dst.Vid = np.copy(self.Vid)
        dst.Iid = np.copy(self.Iid)
        dst.Ae = np.copy(self.Ae)
        dst.vnet = np.copy(self.vnet)
        dst.V_I = np.copy(self.V_I)
        return dst
    def V(self, node):
        return self.vnet[self.nodes_id[node]]
    def I(self, name):
        if name not in self.etype:
            return
        etype = self.etype[name]
        if etype == "R":
            return (self.vnet[self.r[name][0]]-self.vnet[self.r[name][1]])/self.r[name][2]
        if etype == "C":
            return self.V_I[name]
        if etype == "L":
            return self.l[name][3]
        if etype == "V":
            return self.V_I[name]
        if etype == "I":
            return self.isrc[name][2]
        if etype == "E":
            return self.V_I[name]
        if etype == "G":
            return (self.vnet[self.gsrc[name][3]]-self.vnet[self.gsrc[name][4]])*self.gsrc[name][2]
    def attr(self, name):
        if name not in self.etype:
            return
        etype = self.etype[name]
        if etype == "R":
            return self.r[name][2]
        if etype == "C":
            return self.c[name][2]
        if etype == "L":
            return self.l[name][2]
        if etype == "V":
            return self.vsrc[name][2]
        if etype == "I":
            return self.isrc[name][2]
        if etype == "E":
            return self.esrc[name][2]
        if etype == "G":
            return self.gsrc[name][2]
    def ch(self, name, value:float):
        if name not in self.etype:
            return
        etype = self.etype[name]
        if etype == "R":
            src = self.r[name][0]
            dst = self.r[name][1]
            o_value = self.r[name][2]
            self.fG[src, src] -= 1/o_value-1/value
            self.fG[src, dst] += 1/o_value-1/value
            self.fG[dst, src] += 1/o_value-1/value
            self.fG[dst, dst] -= 1/o_value-1/value
            self.sG[self.nmap[src], src] -= 1/o_value-1/value
            self.sG[self.nmap[src], dst] += 1/o_value-1/value
            self.sG[self.nmap[dst], src] += 1/o_value-1/value
            self.sG[self.nmap[dst], dst] -= 1/o_value-1/value
            self.G[self.nmap[src], self.nmap[src]] -= 1/o_value-1/value
            self.G[self.nmap[src], self.nmap[dst]] += 1/o_value-1/value
            self.G[self.nmap[dst], self.nmap[src]] += 1/o_value-1/value
            self.G[self.nmap[dst], self.nmap[dst]] -= 1/o_value-1/value
            self.r[name][2] = value
            return
        if etype == "C":
            src = self.c[name][0]
            dst = self.c[name][1]
            o_value = self.c[name][2]
            oV = self.c[name][3]
            self.c[name][3] *= o_value/value
            dV = self.c[name][3]-oV
            if src in self.vdep[dst]:
                self.Vid[src] += dV
                for x in self.vdep[src]:
                    self.Vid[x] += dV
            else:
                self.Vid[dst] -= dV
                for x in self.vdep[dst]:
                    self.Vid[x] -= dV
            self.c[name][2] = value
            return
        if etype == "L":
            src = self.l[name][0]
            dst = self.l[name][1]
            o_value = self.l[name][2]
            oI = self.l[name][3]
            self.l[name][3] *= o_value/value
            dI = self.l[name][3]-oI
            self.Iid[self.nmap[src]] -= dI
            self.Iid[self.nmap[dst]] += dI
            self.l[name][2] = value
            return
        if etype == "V":
            src = self.vsrc[name][0]
            dst = self.vsrc[name][1]
            dV = value-self.vsrc[name][2]
            if src in self.vdep[dst]:
                self.Vid[src] += dV
                for x in self.vdep[src]:
                    self.Vid[x] += dV
            else:
                self.Vid[dst] -= dV
                for x in self.vdep[dst]:
                    self.Vid[x] -= dV
            self.vsrc[name][2] = value
            return
        if etype == "I":
            src = self.isrc[name][0]
            dst = self.isrc[name][1]
            dI = value-self.isrc[name][3]
            self.Iid[self.nmap[src]] -= dI
            self.Iid[self.nmap[dst]] += dI
            self.isrc[name][2] = value
            return
        if etype == "E":
            src = self.esrc[name][0]
            dst = self.esrc[name][1]
            dvalue = value-self.esrc[name][2]
            np = self.esrc[name][3]
            nn = self.esrc[name][4]
            if dst in self.vdep[src]:
                self.Ve[dst, self.nmap[np]] -= dvalue
                self.Ve[dst, self.nmap[nn]] += dvalue
                for y in self.vdep[dst]:
                    self.Ve[y, self.nmap[np]] -= dvalue
                    self.Ve[y, self.nmap[nn]] += dvalue
            else:
                self.Ve[src, self.nmap[np]] += dvalue
                self.Ve[src, self.nmap[nn]] -= dvalue
                for y in self.vdep[src]:
                    self.Ve[y, self.nmap[np]] += dvalue
                    self.Ve[y, self.nmap[nn]] -= dvalue
            self.esrc[name][2] = value
            self.Ae = np.identity(len(self.nodes))
            for x in self.edep:
                it = self.esrc[x]
                if it[0] in self.vdep[it[1]]:
                    self.Ae[it[0], ::] += it[2]*self.Ae[it[3], ::]
                    self.Ae[it[0], ::] -= it[2]*self.Ae[it[4], ::]
                    for y in self.vdep[it[0]]:
                        self.Ae[y, ::] += it[2]*self.Ae[it[3], ::]
                        self.Ae[y, ::] -= it[2]*self.Ae[it[4], ::]
                else:
                    self.Ae[it[1], ::] -= it[2]*self.Ae[it[3], ::]
                    self.Ae[it[1], ::] += it[2]*self.Ae[it[4], ::]
                    for y in self.vdep[it[1]]:
                        self.Ae[y, ::] -= it[2]*self.Ae[it[3], ::]
                        self.Ae[y, ::] += it[2]*self.Ae[it[4], ::]
            return
        if etype == "G":
            src = self.gsrc[name][0]
            dst = self.gsrc[name][1]
            dvalue = value-self.gsrc[name][2]
            np = self.gsrc[name][3]
            nn = self.gsrc[name][4]
            self.fG[src, np] += dvalue
            self.fG[src, nn] -= dvalue
            self.fG[dst, np] -= dvalue
            self.fG[dst, nn] += dvalue
            self.sG[self.nmap[src], np] += dvalue
            self.sG[self.nmap[src], nn] -= dvalue
            self.sG[self.nmap[dst], np] -= dvalue
            self.sG[self.nmap[dst], nn] += dvalue
            self.G[self.nmap[src], self.nmap[np]] += dvalue
            self.G[self.nmap[src], self.nmap[nn]] -= dvalue
            self.G[self.nmap[dst], self.nmap[np]] -= dvalue
            self.G[self.nmap[dst], self.nmap[nn]] += dvalue
            self.isrc[name][2] = value
            return
class rk4_sp:
    def __init__(self, sp_it):
        self.sp = sp_it
        self.sp_t = sp()
    def run(self, dt, ut, filt, rt=0, rdt=1):
        t = 0
        n = int(ut//dt)
        recode = []
        cnt = 0
        for i in tqdm(range(n)):
            self.sp.run()
            k1C, k1L = self.sp.update_s()
            '''
            self.sp.deepcopy(self.sp_t)
            self.sp_t.update(dt/2, t+dt/2, k1C, k1L)
            self.sp_t.run()
            k2C, k2L = self.sp_t.update_s()
            self.sp.deepcopy(self.sp_t)
            self.sp_t.update(dt/2, t+dt/2, k2C, k2L)
            self.sp_t.run()
            k3C, k3L = self.sp_t.update_s()
            self.sp.deepcopy(self.sp_t)
            self.sp_t.update(dt, t+dt, k3C, k3L)
            self.sp_t.run()
            k4C, k4L = self.sp_t.update_s()
            kC = {}
            kL = {}
            for x in k1C:
                kC[x] = (k1C[x]+2*k2C[x]+2*k3C[x]+k4C[x])/6
            for x in k1L:
                kL[x] = (k1L[x]+2*k2L[x]+2*k3L[x]+k4L[x])/6
            '''
            self.sp.update(dt, t+dt, k1C, k1L)
            t += dt
            if t >= rt:
                cnt += 1
                if cnt == rdt:
                    recode.append(filt(self.sp, t))
                    cnt = 0
        return recode
