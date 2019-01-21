######################################################################
#
#   artemide for python
#           A.Vladimirov (16.01.2019)
#
######################################################################

import artemide

def initialize(order):
    """Initialization of artemide
    
        Argument: order = LO, NLO, NNLO
    """
    artemide.harpy.initialize(order)
    print "Welcome to harpy -- the python interface for artemide"

def setNPparameters(l):
    """Setting NP parameters for the model
                Arguments: (l)
                Argument overloading:
                (integer)        = loads replica
                (array)          = set array on NP parameters
    """
    if isinstance(l,list):
        artemide.harpy.setlambda(l)
    else:
        artemide.harpy.setlambda_byreplica(int(l))

def varyScales(c1,c2,c3,c4):
        """Set new scale variation parameters

                Arguments: (c1,c2,c3,c4)
        """
        artemide.harpy.setscalevariation(c1,c2,c3,c4)

###############################################################################
class TMD5:
        """Static class for TMDPDFs

        TMDs without gluon component
        """
        @staticmethod
        def unpolarizedPDF(x,bt,a=None,b=None,c=None):
                """Unpolarized TMD PDF f_1

                Arguments: (x,b)
                Argument overloading:
                (x,b)            = optimal TMD for hadron=1
                (x,b,h)          = optimal TMD for hadron=h (integer)
                (x,b,mu,zeta)    = TMD at scale (mu,zeta) for hadron=1
                (x,b,mu,zeta,h)  = TMD at scale (mu,zeta) for hadron=h(integer)
                """
                if a==None:                 #(-,?,?)
                    if b==None and c==None: #(-,-,-)
                        return artemide.harpy.utmdpdf_5_optimal(x,bt,1)
                    else:                   #(-,?,?) ?!
                        print "I am not sure how python works. 1"
                        return 0
                elif b==None:               #(a,-,?)
                    if c==None:             #(a,-,-) => a=h
                        return artemide.harpy.utmdpdf_5_optimal(x,bt,int(a))
                    else:                   #(a,-,?) ?!
                        print "I am not sure how python works. 2"
                        return 0
                elif c==None:               #(a,b,-) => (a,b)=(mu,zeta)
                    return artemide.harpy.utmdpdf_5_evolved(x,bt,a,b,1)
                else:
                    return artemide.harpy.utmdpdf_5_evolved(x,bt,a,b,int(c))

###############################################################################
class TMD50:
        """Static class for TMDPDFs

        TMDs with gluon component
        """
        @staticmethod
        def unpolarizedPDF(x,bt,a=None,b=None,c=None):
                """Unpolarized TMD PDF f_1

                Arguments: (x,b)
                Argument overloading:
                (x,b)            = optimal TMD for hadron=1
                (x,b,h)          = optimal TMD for hadron=h (integer)
                (x,b,mu,zeta)    = TMD at scale (mu,zeta) for hadron=1
                (x,b,mu,zeta,h)  = TMD at scale (mu,zeta) for hadron=h(integer)
                """
                if a==None:                 #(-,?,?)
                    if b==None and c==None: #(-,-,-)
                        return artemide.harpy.utmdpdf_50_optimal(x,bt,1)
                    else:                   #(-,?,?) ?!
                        print "I am not sure how python works. 1"
                        return 0
                elif b==None:               #(a,-,?)
                    if c==None:             #(a,-,-) => a=h
                        return artemide.harpy.utmdpdf_50_optimal(x,bt,int(a))
                    else:                   #(a,-,?) ?!
                        print "I am not sure how python works. 2"
                        return 0
                elif c==None:               #(a,b,-) => (a,b)=(mu,zeta)
                    return artemide.harpy.utmdpdf_50_evolved(x,bt,a,b,1)
                else:
                    return artemide.harpy.utmdpdf_50_evolved(x,bt,a,b,int(c))

###############################################################################
class TMD5_inKT:
        """Static class for TMDPDFs

        TMDs without gluon component
        """
        @staticmethod
        def unpolarizedPDF(x,kt,a=None,b=None,c=None):
                """Unpolarized TMD PDF f_1

                Arguments: (x,kt)
                Argument overloading:
                (x,kt)            = optimal TMD for hadron=1
                (x,kt,h)          = optimal TMD for hadron=h (integer)
                (x,kt,mu,zeta)    = TMD at scale (mu,zeta) for hadron=1
                (x,kt,mu,zeta,h)  = TMD at scale (mu,zeta) for hadron=h(integer)
                """
                if a==None:                 #(-,?,?)
                    if b==None and c==None: #(-,-,-)
                        return artemide.harpy.utmdpdf_kt_5_optimal(x,kt,1)
                    else:                   #(-,?,?) ?!
                        print "I am not sure how python works. 1"
                        return 0
                elif b==None:               #(a,-,?)
                    if c==None:             #(a,-,-) => a=h
                        return artemide.harpy.utmdpdf_kt_5_optimal(x,kt,int(a))
                    else:                   #(a,-,?) ?!
                        print "I am not sure how python works. 2"
                        return 0
                elif c==None:               #(a,b,-) => (a,b)=(mu,zeta)
                    return artemide.harpy.utmdpdf_kt_5_evolved(x,kt,a,b,1)
                else:
                    return artemide.harpy.utmdpdf_kt_5_evolved(x,kt,a,b,int(c))

###############################################################################
class TMD50_inKT:
        """Static class for TMDPDFs

        TMDs with gluon component
        """
        @staticmethod
        def unpolarizedPDF(x,kt,a=None,b=None,c=None):
                """Unpolarized TMD PDF f_1

                Arguments: (x,kt)
                Argument overloading:
                (x,kt)            = optimal TMD for hadron=1
                (x,kt,h)          = optimal TMD for hadron=h (integer)
                (x,kt,mu,zeta)    = TMD at scale (mu,zeta) for hadron=1
                (x,kt,mu,zeta,h)  = TMD at scale (mu,zeta) for hadron=h(integer)
                """
                if a==None:                 #(-,?,?)
                    if b==None and c==None: #(-,-,-)
                        return artemide.harpy.utmdpdf_kt_50_optimal(x,kt,1)
                    else:                   #(-,?,?) ?!
                        print "I am not sure how python works. 1"
                        return 0
                elif b==None:               #(a,-,?)
                    if c==None:             #(a,-,-) => a=h
                        return artemide.harpy.utmdpdf_kt_50_optimal(x,kt,int(a))
                    else:                   #(a,-,?) ?!
                        print "I am not sure how python works. 2"
                        return 0
                elif c==None:               #(a,b,-) => (a,b)=(mu,zeta)
                    return artemide.harpy.utmdpdf_kt_50_evolved(x,kt,a,b,1)
                else:
                    return artemide.harpy.utmdpdf_kt_50_evolved(x,kt,a,b,int(c))

###############################################################################
class DY:
        """Static class for evaluation of DY cross-section
        """

        @staticmethod
        def xSec(process,s,qT,Q,y,includeCuts,CutParameters=None,Num=4):
                """Cross-section for DY integrated over bin
                      
                Arguments: (process,s,qT,Q,y,includeCuts,CutParameters=None,Num=4)
                process         = (int, int, int) (see definition in artemide manual)
                s               = Mandelshtan variable s
                qT              = (qT-Min,qT-Max) boundaries of bin in qT
                Q               = (Q-Min,Q-Max) boundaries of bin in Q
                y               = (y-Min,y-Max) boundaries of bin in y
                includeCuts     = True/False to include leptonic cuts
                CutParameters   = (real,real,real,real) must be if includeCuts=True (see definition in artemide manual)
                Num             = even integer, number of section of qt-integration (defaul=4)
                """

                if includeCuts and CutParameters==None:
                    print "ERROR 1: specify the cut parameters"
                    return 0
                
                if includeCuts and (len(CutParameters) is not 4):
                        print "ERROR 2: legnth of CutParameters must 4"
                        return 0
                
                if len(process) is not 3:
                        print "ERROR 3: legnth of process must 3"
                        return 0
                
                if len(qT) is not 2:
                        print "ERROR 4: legnth of qT must 2"
                        return 0
                
                if len(Q) is not 2:
                        print "ERROR 5: legnth of Q must 2"
                        return 0
                        
                if len(y) is not 2:
                        print "ERROR 6: legnth of y must 2"
                        return 0
                
                if not includeCuts:
                        cc=[0,0,0,0]
                else:
                        cc=CutParameters
                
                return artemide.harpy.dy_xsec_single(process,s,qT,Q,y,includeCuts,cc,Num)
        
            
        @staticmethod
        def xSecList(process,s,qT,Q,y,includeCuts,CutParameters):
                return artemide.harpy.dy_xsec_list(process,s,qT,Q,y,includeCuts,CutParameters,len(s))
                                
