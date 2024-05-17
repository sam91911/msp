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
        self.etype = {}
        self.f = []
        self.nodes = ['0']
        self.nodes_id = {'0': 0}
        self.cnt = 0
    def R(self, src, dst, value, name=None, f=None):
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
        if f is not None:
            self.f.append(f)
        return
    def C(self, src, dst, value, ic=0, name=None, f=None):
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
        if f is not None:
            self.f.append(f)
        return
    def L(self, src, dst, value, ic=0, name=None, f=None):
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
        if f is not None:
            self.f.append(f)
        return
    def Vsrc(self, src, dst, value, name=None, f=None):
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
        if f is not None:
            self.f.append(f)
        return
    def Isrc(self, src, dst, value, name=None, f=None):
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
        if f is not None:
            self.f.append(f)
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
        for x in range(len(self.nodes)):
            self.root_node(x)
        for x in self.cnodes:
            if x not in nodemap:
                nodemap[x] = cnt
                cnt += 1
        self.nlen = cnt
        for x in self.cnodes:
            self.nmap.append(nodemap[self.root_node(x)])
        self.vdep = []
        depend = list([i in self.cnodes for i in range(len(self.nodes))])
        while (not all(depend)):
            for x in self.c:
                if (depend[self.c[x][0]])and(not depend[self.c[x][1]]):
                    depend[self.c[x][1]] = True
                    self.vdep.append(("C", x, True))
                if (depend[self.c[x][1]])and(not depend[self.c[x][0]]):
                    depend[self.c[x][0]] = True
                    self.vdep.append(("C", x, False))
            for x in self.vsrc:
                if (depend[self.vsrc[x][0]])and(not depend[self.vsrc[x][1]]):
                    depend[self.vsrc[x][1]] = True
                    self.vdep.append(("V", x, True))
                if (depend[self.vsrc[x][1]])and(not depend[self.vsrc[x][0]]):
                    depend[self.vsrc[x][0]] = True
                    self.vdep.append(("V", x, False))
        self.videp = []
        depend = [0]*len(self.nodes)
        for x in self.c.values():
            depend[x[0]] += 1
            depend[x[1]] += 1
        for x in self.vsrc.values():
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
        print("end setup")
        return
    def get_vnet(self):
        vnet = [0]*len(self.nodes)
        for x in self.vdep:
            name = x[1]
            if x[0] == "C":
                it = self.c[name]
                if x[2]:
                    vnet[it[1]] = vnet[it[0]] - it[3]
                else:
                    vnet[it[0]] = vnet[it[1]] + it[3]
            if x[0] == "V":
                it = self.vsrc[name]
                if x[2]:
                    vnet[it[1]] = vnet[it[0]] - it[2]
                else:
                    vnet[it[0]] = vnet[it[1]] + it[2]
        return vnet
    def get_fG(self):
        fG = np.zeros((len(self.nodes), len(self.nodes)))
        for x in self.r.values():
            fG[x[0], x[0]] += 1/x[2]
            fG[x[0], x[1]] -= 1/x[2]
            fG[x[1], x[1]] += 1/x[2]
            fG[x[1], x[0]] -= 1/x[2]
        return fG
    def get_sG(self):
        sG = np.zeros((self.nlen, len(self.nodes)))
        for x in self.r.values():
            sG[self.nmap[x[0]], x[0]] += 1/x[2]
            sG[self.nmap[x[0]], x[1]] -= 1/x[2]
            sG[self.nmap[x[1]], x[1]] += 1/x[2]
            sG[self.nmap[x[1]], x[0]] -= 1/x[2]
        return sG
    def get_G(self):
        G = np.zeros((self.nlen, self.nlen))
        for x in self.r.values():
            G[self.nmap[x[0]], self.nmap[x[0]]] += 1/x[2]
            G[self.nmap[x[0]], self.nmap[x[1]]] -= 1/x[2]
            G[self.nmap[x[1]], self.nmap[x[1]]] += 1/x[2]
            G[self.nmap[x[1]], self.nmap[x[0]]] -= 1/x[2]
        return G
    def get_Id(self):
        I = np.zeros((self.nlen,))
        for x in self.l.values():
            I[self.nmap[x[0]]] -= x[3]
            I[self.nmap[x[1]]] += x[3]
        for x in self.isrc.values():
            I[self.nmap[x[0]]] -= x[2]
            I[self.nmap[x[1]]] += x[2]
        vnet = self.get_vnet()
        sG = self.get_sG()
        for i in range(len(vnet)):
            I -= sG[::, i]*vnet[i]
        return I
    def get_V(self):
        G = self.get_G()[::, 1::]
        I = self.get_Id()
        V = np.linalg.lstsq(G, I, rcond = None)[0]
        V = [0]+V.tolist()
        vnet = self.get_vnet()
        for i in range(len(vnet)):
            vnet[i] += V[self.nmap[i]]
        return vnet
    def get_I(self, V = None):
        I = np.zeros((len(self.nodes),))
        fG = self.get_fG()
        for x in self.l.values():
            I[x[0]] -= x[3]
            I[x[1]] += x[3]
        for x in self.isrc.values():
            I[x[0]] -= x[2]
            I[x[1]] += x[2]
        if V is None:
            V = self.get_V()
        for i in range(len(V)):
            I -= fG[::, i]*V[i]
        C_I = {}
        V_I = {}
        for x in self.videp:
            if x[0] == "C":
                src = self.c[x[1]][0]
                dst = self.c[x[1]][1]
                if x[2]:
                    ei = I[src]
                else:
                    ei = -I[dst]
                C_I[x[1]] = ei
                I[src] -= ei
                I[dst] += ei
            if x[0] == "V":
                src = self.vsrc[x[1]][0]
                dst = self.vsrc[x[1]][1]
                if x[2]:
                    ei = I[src]
                else:
                    ei = -I[dst]
                V_I[x[1]] = -ei
                I[src] -= ei
                I[dst] += ei
        return C_I, V_I
    def update_s(self):
        V = self.get_V()
        C_I, _ = self.get_I(V)
        C_dv = {}
        L_di = {}
        for x in self.c:
            C_dv[x] = C_I[x]/self.c[x][2]
        for x in self.l:
            dv = V[self.l[x][0]] - V[self.l[x][1]]
            L_di[x] = dv/self.l[x][2]
        return C_dv, L_di
    def update(self, dt, t, C_dv, L_di):
        V = self.get_V()
        C_I, V_I = self.get_I(V)
        for x in self.c:
            self.c[x][3] += C_dv[x]*dt
        for x in self.l:
            self.l[x][3] += L_di[x]*dt
        for f in self.f:
            f(self, (V, C_I, V_I), dt, t)
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
        dst.cnodes = copy.deepcopy(self.cnodes)
        dst.rank = copy.deepcopy(self.rank)
        dst.nmap = copy.deepcopy(self.nmap)
        dst.nlen = copy.deepcopy(self.nlen)
        dst.cnt = copy.deepcopy(self.cnt)
        return dst
    def V(self, node, state):
        return state[0][self.nodes_id[node]]
    def I(self, name, state):
        if name not in self.etype:
            return
        etype = self.etype[name]
        if etype == "R":
            return (state[0][self.R[name][0]]-state[0][self.R[name][1]])/self.R[name][2]
        if etype == "C":
            return state[1][name]
        if etype == "L":
            return self.l[name][3]
        if etype == "V":
            return state[2][name]
        if etype == "I":
            return self.isrc[name][2]
    def ch(self, name, value:float):
        if name not in self.etype:
            return
        etype = self.etype[name]
        if etype == "R":
            self.r[name][2] = value
            return
        if etype == "C":
            self.c[name][2] = value
            return
        if etype == "L":
            self.l[name][2] = value
            return
        if etype == "V":
            self.vsrc[name][2] = value
            return
        if etype == "I":
            self.isrc[name][2] = value
            return
class rk4_sp:
    def __init__(self, sp_it):
        self.sp = sp_it
        self.sp_t = self.sp.deepcopy()
    def run(self, dt, ut, rt, filt):
        t = 0
        n = int(ut//dt)
        recode = []
        for i in tqdm(range(n)):
            k1C, k1L = self.sp.update_s()
            self.sp.deepcopy(self.sp_t)
            self.sp_t.update(dt/2, t+dt/2, k1C, k1L)
            k2C, k2L = self.sp_t.update_s()
            self.sp.deepcopy(self.sp_t)
            self.sp_t.update(dt/2, t+dt/2, k2C, k2L)
            k3C, k3L = self.sp_t.update_s()
            self.sp.deepcopy(self.sp_t)
            self.sp_t.update(dt, t+dt, k3C, k3L)
            k4C, k4L = self.sp_t.update_s()
            kC = {}
            kL = {}
            for x in k1C:
                kC[x] = (k1C[x]+2*k2C[x]+2*k3C[x]+k4C[x])/6
            for x in k1L:
                kL[x] = (k1L[x]+2*k2L[x]+2*k3L[x]+k4L[x])/6
            self.sp.update(dt, t+dt, kC, kL)
            t += dt
            if t >= rt:
                V = self.sp.get_V()
                C_I, V_I = self.sp.get_I(V)
                recode.append(filt(self.sp, (V, C_I, V_I), t))
        return recode
