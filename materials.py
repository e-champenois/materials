import sympy  as sp
import numpy as np
from scipy.constants import c

l = sp.symbols('l')
w = sp.symbols('w')
th = sp.symbols('th')

class Material():
    '''
    n^2 - 1 = n0 + l^2 * A / (l^2 - B^2) + C / (l^2 - D^2) + E * l2
    Data from refractiveindex.info
    '''
    def __init__(self, A=np.array([]), B=np.array([]),
                       C=np.array([]), D=np.array([]),
                       E=np.array([]), n0=0):
        
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E
        self.n0 = n0
        
        self.nfoo = 1 + n0
        for Ai, Bi in zip(self.A, self.B):
            self.nfoo += l**2 * Ai / (l**2 - Bi**2)
        for Ci, Di in zip(self.C, self.D):
            self.nfoo += Ci / (l**2 - Di**2)
        for Ei in self.E:
            self.nfoo += Ei * l**2
        self.nfoo = sp.sqrt(self.nfoo)
        
    def n(self, ls):
        self.n = sp.lambdify(l, self.nfoo, np)
        return self.n(ls)
        
    def ng(self, ls):
        self.make_dkfoo()
        self.ng = sp.lambdify(l, self.dkfoo.subs(w, 2*np.pi*c/l) * c, np)
        return self.ng(ls)
        
    def GVD(self, ls):
        self.make_GVDfoo()
        self.GVD = sp.lambdify(l, self.GVDfoo.subs(w, 2*np.pi*c/l), np)
        return self.GVD(ls)
            
    def TOD(self, ls):
        self.make_TODfoo()
        self.TOD = sp.lambdify(l, self.TODfoo.subs(w, 2*np.pi*c/l), np)
        return self.TOD(ls)
    
    def HOD(self, ls, order=4):
        self.make_kfoo()
        hod = self.kfoo
        for i in range(order):
            hod = hod.diff(w)
        return sp.lambdify(l, hod.subs(w, 2*np.pi*c/l), np)(ls)
        
    def make_kfoo(self):
        if not hasattr(self, 'kfoo'):
            self.kfoo = w/c * self.nfoo.subs(l, 2*np.pi*c/w)
                    
    def make_dkfoo(self):
        if not hasattr(self, 'dkfoo'):
            self.make_kfoo()
            self.dkfoo = self.kfoo.diff(w)
                
    def make_GVDfoo(self):
        if not hasattr(self, 'GVDfoo'):
            self.make_dkfoo()
            self.GVDfoo = self.dkfoo.diff(w)
                
    def make_TODfoo(self):
        if not hasattr(self, 'TODfoo'):
            self.make_GVDfoo()
            self.TODfoo = self.GVDfoo.diff(w)
            
class Birefringent():
    '''
    nth = 1 / sqrt( ( sin(th) / ne )**2 + ( cos(th) / no )**2 )
    '''
    def __init__(self, Mato, Mate):
        
        self.Mato = Mato
        self.Mate = Mate
        self.nfooth = 1 / sp.sqrt( (sp.sin(th)/Mate.nfoo)**2 + (sp.cos(th)/Mato.nfoo)**2 )
        
        self.no = self.Mato.n
        self.ne = self.Mate.n
        
        self.ngo = self.Mato.ng
        self.nge = self.Mate.ng
        
        self.GVDo = self.Mato.GVD
        self.GVDe = self.Mate.GVD
        
        self.TODo = self.Mato.TOD
        self.TODe = self.Mate.TOD
        
        self.HODo = self.Mato.HOD
        self.HODe = self.Mate.HOD
        
    def nth(self, ls, ths):
        self.nth = sp.lambdify([l, th], self.nfooth, np)
        return self.nth(ls, ths)

    def ngth(self, ls, ths):
        self.make_dkfooth()
        self.ngth = sp.lambdify([l, th], self.dkfooth.subs(w, 2*np.pi*c/l) * c, np)
        return self.ngth(ls, ths)

    def GVDth(self, ls, ths):
        self.make_GVDfooth()
        self.GVD = sp.lambdify([l, th], self.GVDfooth.subs(w, 2*np.pi*c/l), np)
        return self.GVD(ls, ths)

    def TODth(self, ls, ths):
        self.make_TODfooth()
        self.TOD = sp.lambdify([l, th], self.TODfooth.subs(w, 2*np.pi*c/l), np)
        return self.TOD(ls, ths)

    def make_kfooth(self):
        if not hasattr(self, 'kfooth'):
            self.kfooth = w/c * self.nfooth.subs(l, 2*np.pi*c/w)

    def make_dkfooth(self):
        if not hasattr(self, 'dkfooth'):
            self.make_kfooth()
            self.dkfooth = self.kfooth.diff(w)

    def make_GVDfooth(self):
        if not hasattr(self, 'GVDfooth'):
            self.make_dkfooth()
            self.GVDfooth = self.dkfooth.diff(w)

    def make_TODfooth(self):
        if not hasattr(self, 'TODfooth'):
            self.make_GVDfooth()
            self.TODfooth = self.GVDfooth.diff(w)
        
C = np.array([0.01878]) * 1e-12
D = np.sqrt(np.array([0.01822])) * 1e-6
E = np.array([-0.01354]) * 1e12
n0 = 2.7471 - 1
aBBOo = Material(C=C, D=D, E=E, n0=n0)

C = np.array([0.01224]) * 1e-12
D = np.sqrt(np.array([0.01667])) * 1e-6
E = np.array([-0.01516]) * 1e12
n0 = 2.3715 - 1
aBBOe = Material(C=C, D=D, E=E, n0=n0)

aBBO = Birefringent(aBBOo, aBBOe)

A = np.array([15.102464])
B = np.sqrt(np.array([400])) * 1e-6
C = np.array([0.011125165]) * 1e-12
D = np.sqrt(np.array([0.01325366])) * 1e-6
n0 = 2.302842 - 1
ADPo = Material(A=A, B=B, C=C, D=D, n0=n0)

A = np.array([5.919896])
B = np.sqrt(np.array([400])) * 1e-6
C = np.array([0.009616676]) * 1e-12
D = np.sqrt(np.array([0.01298912])) * 1e-6
n0 = 2.163510 - 1
ADPe = Material(A=A, B=B, C=C, D=D, n0=n0)

ADP = Birefringent(ADPo, ADPe)

C = np.array([0.2311]) * 1e-12
D = np.sqrt(np.array([0.0688])) * 1e-6
E = np.array([-0.00257]) * 1e12
n0 = 5.7975 - 1
AgGaS2o = Material(C=C, D=D, E=E, n0=n0)

C = np.array([0.2230]) * 1e-12
D = np.sqrt(np.array([0.0946])) * 1e-6
E = np.array([-0.00261]) * 1e12
n0 = 5.5436 - 1
AgGaS2e = Material(C=C, D=D, E=E, n0=n0)

AgGaS2 = Birefringent(AgGaS2o, AgGaS2e)

A = np.array([0.8107, 0.19652, 4.52469])
B = np.array([0.10065, 29.87, 53.82]) * 1e-6
n0 = 0.33973
BaF2 = Material(A=A, B=B, n0=n0)

A = np.array([0.90291, 0.83155, 0.76536])
B = np.sqrt(np.array([0.003926, 0.018786, 60.01])) * 1e-6
BBOo = Material(A=A, B=B)

A = np.array([1.151075, 0.21803, 0.656])
B = np.sqrt(np.array([0.007142, 0.02259, 263])) * 1e-6
BBOe = Material(A=A, B=B)

BBO = Birefringent(BBOo, BBOe)

A = np.array([1.03961212, 0.231792344, 1.01046945])
B = np.sqrt(np.array([0.00600069867, 0.0200179144, 103.560653])) * 1e-6
BK7 = Material(A=A, B=B)

A = np.array([0.96464345, 1.82831454])
B = np.sqrt(np.array([1.94325203e-2, 120])) * 1e-6
n0 = 0.73358749
CaCO3o = Material(A=A, B=B, n0=n0)

A = np.array([0.82427830, 0.14429128])
B = np.sqrt(np.array([1.06689543e-2, 120])) * 1e-6
n0 = 0.35859695
CaCO3e = Material(A=A, B=B, n0=n0)

CaCO3 = Birefringent(CaCO3o, CaCO3e)

A = np.array([0.5675888, 0.4710914, 3.8484723])
B = np.array([0.050263605, 0.1003909, 34.649040]) * 1e-6
CaF2 = Material(A=A, B=B)

A = np.array([0.6961663, 0.4079426, 0.8974794])
B = np.array([0.0684043, 0.1162414, 9.896161]) * 1e-6
FS = Material(A=A, B=B)

A = np.array([0.683740494,0.420323613,0.58502748])
B = np.sqrt(np.array([0.00460352869,0.0133968856,64.4932732])) * 1e-6
FSir = Material(A=A, B=B)

# https://www.heraeus.com/media/media/hca/doc_hca/products_and_solutions_8/optics/Data_and_Properties_Optics_fused_silica_EN.pdf
A = np.array([4.76523070e-1, 6.27786368e-1, 8.72274404e-1])
B = np.sqrt(np.array([2.84888095e-3, 1.18369052e-2, 9.56856012e1])) * 1e-6
Infrasil = Material(A=A, B=B)

A = np.array([0.79221, 0.01981, 0.15587, 0.17673, 2.06217])
B = np.array([0.146, 0.173, 0.187, 60.61, 87.72]) * 1e-6
n0 = 0.39408
KBr = Material(A=A, B=B, n0=n0)

A = np.array([13.00522])
B = np.sqrt(np.array([400])) * 1e-6
C = np.array([0.01008956]) * 1e-12
D = np.sqrt(np.array([0.0129426])) * 1e-6
n0 = 2.259276 - 1
KDPo = Material(A=A, B=B, C=C, D=D, n0=n0)

A = np.array([3.2279924])
B = np.sqrt(np.array([400])) * 1e-6
C = np.array([0.008637494]) * 1e-12
D = np.sqrt(np.array([0.0122810])) * 1e-6
n0 = 2.132668 - 1
KDPe = Material(A=A, B=B, C=C, D=D, n0=n0)

KDP = Birefringent(KDPo, KDPe)

A = np.array([1.8293958, 1.6675593, 1.1210424, 0.04513366, 12.380234])
B = np.sqrt(np.array([0.0225, 0.0625, 0.1225, 0.2025, 27089.737])) * 1e-6
KRS5 = Material(A=A, B=B)

C = np.array([0.04140, 9.35522]) * 1e-12
D = np.sqrt(np.array([0.03978, 31.45571])) * 1e-6
n0 = 3.29100 - 1
KTPa = Material(C=C, D=D, n0=n0)

C = np.array([0.04341, 16.98825]) * 1e-12
D = np.sqrt(np.array([0.04597, 39.43799])) * 1e-6
n0 = 3.45018 - 1
KTPb = Material(C=C, D=D, n0=n0)

C = np.array([0.06206, 110.80672]) * 1e-12
D = np.sqrt(np.array([0.04763, 86.12171])) * 1e-6
n0 = 4.59423 - 1
KTPg = Material(C=C, D=D, n0=n0)

A = np.array([1.37623, 1.06745])
B = np.sqrt(np.array([0.0350832, 169])) * 1e-6
n0 = 1.03132
LiIO3o = Material(A=A, B=B, n0=n0)

A = np.array([1.08807, 0.554582])
B = np.sqrt(np.array([0.0313810, 158.76])) * 1e-6
n0 = 0.83086
LiIO3e = Material(A=A, B=B, n0=n0)

LiIO3 = Birefringent(LiIO3o, LiIO3e)

A = np.array([2.6734, 1.2290, 12.614])
B = np.sqrt(np.array([0.01764, 0.05914, 474.6])) * 1e-6
LiNbO3o = Material(A=A, B=B)

A = np.array([2.9804, 0.5981, 8.9543])
B = np.sqrt(np.array([0.02047, 0.0666, 416.08])) * 1e-6
LiNbO3e = Material(A=A, B=B)

LiNbO3 = Birefringent(LiNbO3o, LiNbO3e)

A = np.array([0.48755108, 0.39875031, 2.3120353])
B = np.array([0.04338408, 0.09461442, 23.793604]) * 1e-6
MgF2o = Material(A=A, B=B)

A = np.array([0.41344023, 0.50497499, 2.4904862])
B = np.array([0.03684262, 0.09076162, 23.771995]) * 1e-6
MgF2e = Material(A=A, B=B)

MgF2 = Birefringent(MgF2o, MgF2e)

A = np.array([1.4313493, 0.65054713, 5.3414021])
B = np.array([0.0726631, 0.1193242, 18.028251]) * 1e-6
Sappho = Material(A=A, B=B)

A = np.array([1.5039759, 0.55069141, 6.5927379])
B = np.array([0.0740288, 0.1216529, 20.072248]) * 1e-6
Sapphe = Material(A=A, B=B)

Sapph = Birefringent(Sappho, Sapphe)

A = np.array([1.55912932, 0.284246288, 0.968842926])
B = np.sqrt(np.array([0.0121481001, 0.0534549042, 112.174809])) * 1e-6
SF1 = Material(A=A, B=B)

A = np.array([1.40301821, 0.231767504, 0.939056586])
B = np.sqrt(np.array([0.0105795466, 0.0493226978, 112.405955])) * 1e-6
SF2 = Material(A=A, B=B)

A = np.array([1.61957826, 0.339493189, 1.02566931])
B = np.sqrt(np.array([0.0125502104, 0.0544559822, 117.652222])) * 1e-6
SF4 = Material(A=A, B=B)

A = np.array([1.46141885, 0.247713019, 0.949995832])
B = np.sqrt(np.array([0.0111826126, 0.0508594669, 112.041888])) * 1e-6
SF5 = Material(A=A, B=B)

A = np.array([1.72448482, 0.390104889, 1.04572858])
B = np.sqrt(np.array([0.0134871947, 0.0569318095, 118.557185])) * 1e-6
SF6 = Material(A=A, B=B)

A = np.array([1.61625977, 0.259229334, 1.07762317])
B = np.sqrt(np.array([0.0127534559, 0.0581983954, 116.60768])) * 1e-6
SF10 = Material(A=A, B=B)

A = np.array([1.73848403, 0.311168974, 1.17490871])
B = np.sqrt(np.array([0.0136068604, 0.0615960463, 121.922711])) * 1e-6
SF11 = Material(A=A, B=B)

A = np.array([1.69182538, 0.285919934, 1.12595145])
B = np.sqrt(np.array([0.0133151542, 0.0612647445, 118.405242])) * 1e-6
SF14 = Material(A=A, B=B)

A = np.array([1.53925927, 0.247620926, 1.03816409])
B = np.sqrt(np.array([0.0119307961, 0.0556077536, 116.416747])) * 1e-6
SF15 = Material(A=A, B=B)

A = np.array([1.81651371, 0.428893641, 1.07186278])
B = np.sqrt(np.array([0.0143704198, 0.0592801172, 121.419942])) * 1e-6
SF57 = Material(A=A, B=B)

A = np.array([2.07842233, 0.407120032, 1.76711292])
B = np.sqrt(np.array([0.0180875134, 0.0679493572, 215.266127])) * 1e-6
SF66 = Material(A=A, B=B)

C = np.array([1, 0.004482633]) * 1e-12
D = np.array([0, 1.108205]) * 1e-6
n0 = 11.67316 - 1
Si = Material(C=C, D=D, n0=n0)

A = np.array([1.07044083, 1.10202242])
B = np.sqrt(np.array([1.00585997e-2, 100])) * 1e-6
n0 = 0.28604141
SiO2o = Material(A=A, B=B, n0=n0)

A = np.array([1.09509924, 1.15662475])
B = np.sqrt(np.array([1.02101864e-2, 100])) * 1e-6
n0 = 0.28851804
SiO2e = Material(A=A, B=B, n0=n0)

SiO2 = Birefringent(SiO2o, SiO2e)

A = np.array([0.7097, 0.1788, 3.8796])
B = np.array([0.09597, 26.03, 45.60]) * 1e-6
n0 = 0.33973
SrF2 = Material(A=A, B=B, n0=n0)

C = np.array([0.2441]) * 1e-12
D = np.sqrt(np.array([0.0803])) * 1e-6
n0 = 5.912 - 1
TiO2o = Material(C=C, D=D, n0=n0)

C = np.array([0.3322]) * 1e-12
D = np.sqrt(np.array([0.0843])) * 1e-6
n0 = 7.197 - 1
TiO2e = Material(C=C, D=D, n0=n0)

TiO2 = Birefringent(TiO2o, TiO2e)

A = np.array([2.282, 3.27644])
B = np.sqrt(np.array([0.01185, 282.734])) * 1e-6
YAG = Material(A=A, B=B)

A = np.array([4.45813734, 0.467216334, 2.89566290])
B = np.array([0.200859853, 0.391371166, 47.1362108]) * 1e-6
ZnSe = Material(A=A, B=B)

materials = {'aBBOo':aBBOo, 'aBBOe':aBBOe, 'aBBO':aBBO, 'ADPo':ADPo, 'ADPe':ADPe, 'ADP':ADP,
             'AgGaS2o':AgGaS2o, 'AgGaS2e':AgGaS2e, 'AgGaS2':AgGaS2, 'BaF2':BaF2,
             'BBOo':BBOo, 'BBOe':BBOe, 'BBO':BBO, 'BK7':BK7,
             'CaCO3o':CaCO3o, 'CaCO3e':CaCO3e, 'CaCO3':CaCO3, 'CaF2':CaF2, 'FS':FS,
             'FSir':FSir, 'Infrasil':Infrasil, 'KBr':KBr, 'KDPo':KDPo, 'KDPe':KDPe, 'KDP':KDP,
             'KRS5':KRS5, 'KTPa':KTPa, 'KTPb':KTPb, 'KTPg':KTPg, 'LiIO3o':LiIO3o, 'LiIO3e':LiIO3e, 'LiIO3':LiIO3,
             'LiNbO3o':LiNbO3o, 'LiNbO3e':LiNbO3e, 'LiNbO3':LiNbO3, 'MgF2o':MgF2o, 'MgF2e':MgF2e, 'MgF2':MgF2,
             'Sappho':Sappho, 'Sapphe':Sapphe, 'Sapph':Sapph,
             'SF1':SF1, 'SF2':SF2, 'SF4':SF4, 'SF5':SF5, 'SF6':SF6,
             'SF10':SF10, 'SF11':SF11, 'SF14':SF14, 'SF15':SF15,
             'SF57':SF57, 'SF66':SF66, 'Si':Si, 'SiO2o':SiO2o, 'SiO2e':SiO2e, 'SiO2':SiO2,
             'SrF2':SrF2, 'TiO2o':TiO2o, 'TiO2e':TiO2e, 'TiO2':TiO2, 'YAG':YAG, 'ZnSe':ZnSe}
