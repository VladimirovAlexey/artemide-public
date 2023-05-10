######################################################################
#
#   artemide for python2
#           A.Vladimirov (16.01.2019)
#
######################################################################

import artemide
import numpy

def initialize(fileName):
    """Initialization of artemide
    
        Argument: order = LO, NLO, NNLO
    """
    artemide.harpy.initialize(fileName)
    if artemide.harpy.started:
        pass
    else:
        print "Welcome to harpy -- the python interface for artemide"

def ShowStatistics():
    artemide.harpy.showstatistics()

def setNPparameters(l):
    """Setting NP parameters for the model
                Arguments: (l)
                Argument overloading:
                (integer)        = loads replica
                (array)          = set array on NP parameters
    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_main(numpy.asfortranarray(l))
    else:
        print 'ERROR: argument must be list'
        
def setNPparameters_TMDR(l):
    """Setting NP parameters for the model of TMDR
                Arguments: (l)
                Argument overloading:
                (integer)        = loads replica
                (array)          = set array on NP parameters
    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_tmdr(numpy.asfortranarray(l))
    else:
        artemide.harpy.setreplica_tmdr(int(l))

def setNPparameters_uTMDPDF(l):
    """Setting NP parameters for the model of uTMDPDF
                Arguments: (l)
                Argument overloading:
                (integer)        = loads replica
                (array)          = set array on NP parameters
    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_utmdpdf(numpy.asfortranarray(l))
    else:
        artemide.harpy.setreplica_utmdpdf(int(l))
        
        
def setNPparameters_uTMDFF(l):
    """Setting NP parameters for the model of uTMDFF
                Arguments: (l)
                Argument overloading:
                (integer)        = loads replica
                (array)          = set array on NP parameters
    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_utmdff(numpy.asfortranarray(l))
    else:
        artemide.harpy.setreplica_utmdff(int(l))
        
def setNPparameters_lpTMDPDF(l):
    """Setting NP parameters for the model of lpTMDPDF
                Arguments: (l)
                Argument overloading:
                (integer)        = loads replica
                (array)          = set array on NP parameters
    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_lptmdpdf(numpy.asfortranarray(l))
    else:
        artemide.harpy.setreplica_lptmdpdf(int(l))

def setNPparameters_SiversTMDPDF(l):
    """Setting NP parameters for the model of SiversTMDPDF
                Arguments: (l)
                Argument overloading:
                (integer)        = loads replica
                (array)          = set array on NP parameters
    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_siverstmdpdf(numpy.asfortranarray(l))
    else:
        artemide.harpy.setreplica_siverstmdpdf(int(l))
        
def setNPparameters_wgtTMDPDF(l):
    """Setting NP parameters for the model of wgtTMDPDF
                Arguments: (l)
                Argument overloading:
                (integer)        = loads replica
                (array)          = set array on NP parameters
    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_wgttmdpdf(numpy.asfortranarray(l))
    else:
        artemide.harpy.setreplica_wgttmdpdf(int(l))


def varyScales(c1,c2,c3,c4):
        """Set new scale variation parameters

                Arguments: (c1,c2,c3,c4)
        """
        artemide.harpy.setscalevariation(c1,c2,c3,c4)

def _IsKinematicProper(s,qT,Q,y):
    """ Checks the point for the proper kinematics
    Especially for the correct domain of X.
    """
    gridX=0.00001
    if qT[0]>qT[1]:
        print 'Wrong order of qT'
        return False
    if Q[0]>Q[1]:
        print 'Wrong order of Q'
        return False
    if y[0]>y[1]:
        print 'Wrong order of y'
        return False
    if qT[1]>Q[0]:
        print 'qT (',qT[1],') > Q(', Q[0],')'
        return False
    if Q[1]>s:
        print 'Q (',Q[1],') > s(', s,')'
        return False
    
    x1x2=(Q[0]**2+qT[0]**2)/s
    ymax=-numpy.log(numpy.sqrt(x1x2))
    ymin=-ymax
    
    if y[1]<ymin or y[0]>ymax:
        print 'y ',y, 'is outside physical region ',[ymin,ymax]
        return False
    
    if y[1]>ymax:
        if x1x2<gridX:
            print 'x outside of the grid'
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(-y[1])<gridX:
            print 'x outside of the grid'
            return False
    
    if y[0]<ymin:
        if x1x2<gridX:
            print 'x outside of the grid'
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(y[0])<gridX:
            print 'x outside of the grid'
            return False
     
    x1x2=(Q[0]**2+qT[0]**2)/s
    ymax=-numpy.log(numpy.sqrt(x1x2))
    ymin=-ymax
    
    if y[1]<ymin or y[0]>ymax:
        print 'y ',y, 'is outside physical region ',[ymin,ymax]
        return False
    
    if y[1]>ymax:
        if x1x2<gridX:
            print 'x outside of the grid'
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(-y[1])<gridX:
            print 'x outside of the grid'
            return False
    
    if y[0]<ymin:
        if x1x2<gridX:
            print 'x outside of the grid'
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(y[0])<gridX:
            print 'x outside of the grid'
            return False
    
    return True

def setPDFreplica(n):
    """Changes the replica for PDF input.
    
        This is a temporary function will be changed in future versions
    """
    artemide.harpy.setpdfreplica(n)
    

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
        def xSec(process,s,qT,Q,y,includeCuts,CutParameters=None,Num=None):
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

#                if includeCuts and CutParameters==None:
#                    print "ERROR 1: specify the cut parameters"
#                    return 0
#                
#                if includeCuts and (len(CutParameters) != 4):
#                        print "ERROR 2: legnth of CutParameters must 4"
#                        return 0
#                
#                if len(process) != 3:
#                        print "ERROR 3: legnth of process must 3"
#                        return 0
#                
#                if len(qT) != 2:
#                        print "ERROR 4: legnth of qT must 2"
#                        return 0
#                
#                if len(Q) != 2:
#                        print "ERROR 5: legnth of Q must 2"
#                        return 0
#                        
#                if len(y) != 2:
#                        print "ERROR 6: legnth of y must 2"
#                        return 0
                
                if not includeCuts:
                        cc=[0,0,0,0]
                else:
                        cc=CutParameters
                
                if Num==None:
                    return artemide.harpy.dy_xsec_single(\
                        numpy.asfortranarray(process),\
                        s,\
                        numpy.asfortranarray(qT),\
                        numpy.asfortranarray(Q),\
                        numpy.asfortranarray(y),\
                        includeCuts,\
                        numpy.asfortranarray(cc))
                else:
                    return artemide.harpy.dy_xsec_single(\
                        numpy.asfortranarray(process),\
                        s,\
                        numpy.asfortranarray(qT),\
                        numpy.asfortranarray(Q),\
                        numpy.asfortranarray(y),\
                        includeCuts,\
                        numpy.asfortranarray(cc),\
                        Num)
                    
            
        @staticmethod
        def xSecList(process,s,qT,Q,y,includeCuts,CutParameters):
            return artemide.harpy.dy_xsec_list(numpy.asfortranarray(process),\
                                               numpy.asfortranarray(s),\
                                               numpy.asfortranarray(qT),\
                                               numpy.asfortranarray(Q),\
                                               numpy.asfortranarray(y),\
                                               numpy.asfortranarray(includeCuts),\
                                               numpy.asfortranarray(CutParameters),\
                                               len(s))
        @staticmethod
        def xSecListBINLESS(process,s,qT,Q,y,includeCuts,CutParameters):
            """ The evaluation of cross-section at a single point. 
            """
            #print len(s)
            
            return artemide.harpy.dy_xsec_binless_list(numpy.asfortranarray(process),\
                                               numpy.asfortranarray(s),\
                                               numpy.asfortranarray(qT),\
                                               numpy.asfortranarray(Q),\
                                               numpy.asfortranarray(y),\
                                               numpy.asfortranarray(includeCuts),\
                                               numpy.asfortranarray(CutParameters),\
                                               len(s))
                                
###############################################################################
class SIDIS:
        """Static class for evaluation of SIDIS cross-section
        """

        @staticmethod
        def xSec(process,s,pT,z,x,Q,includeCuts,CutParameters=None,masses=[0.938,0.130]):
                """Cross-section for DY integrated over bin
                      
                Arguments: (process,s,qT,z,x,Q,includeCuts,CutParameters=None)
                process         = (int, int, int) (see definition in artemide manual)
                s               = Mandelshtan variable s
                pT              = (pT-Min,pT-Max) boundaries of bin in pT
                z               = (z-Min,z-Max) boundaries of bin in z
                x               = (x-Min,x-Max) boundaries of bin in x
                Q               = (Q-Min,Q-Max) boundaries of bin in Q                
                includeCuts     = True/False to include leptonic cuts
                CutParameters   = (real,real,real,real) must be if includeCuts=True (see definition in artemide manual)
                """
                
                if not includeCuts:
                        cc=[0,0,0,100]
                else:
                        cc=CutParameters
                
                
                return artemide.harpy.sidis_xsec_single_withmasses(\
                    numpy.asfortranarray(process),\
                    s,\
                    numpy.asfortranarray(pT),\
                    numpy.asfortranarray(z),\
                    numpy.asfortranarray(x),\
                    numpy.asfortranarray(Q),\
                    includeCuts,\
                    numpy.asfortranarray(cc),
                    numpy.asfortranarray(masses))
                    
            
        @staticmethod
        def xSecList(process,s,pT,z,x,Q,includeCuts,CutParameters,masses):  
            return artemide.harpy.sidis_xsec_list(numpy.asfortranarray(process),\
                                               numpy.asfortranarray(s),\
                                               numpy.asfortranarray(pT),\
                                               numpy.asfortranarray(z),\
                                               numpy.asfortranarray(x),\
                                               numpy.asfortranarray(Q),\
                                               numpy.asfortranarray(includeCuts),\
                                               numpy.asfortranarray(CutParameters),\
                                               numpy.asfortranarray(masses),\
                                               len(s))
            
        @staticmethod
        def xSecListBINLESS(process,s,pT,z,x,Q,masses):  
            """ The evaluation of cross-section at a single point. 
            
            No binning effects.
            Consiquetly, no cuts.
            """
            return artemide.harpy.sidis_xsec_binless_list(numpy.asfortranarray(process),\
                                               numpy.asfortranarray(s),\
                                               numpy.asfortranarray(pT),\
                                               numpy.asfortranarray(z),\
                                               numpy.asfortranarray(x),\
                                               numpy.asfortranarray(Q),\
                                               numpy.asfortranarray(masses),\
                                               len(s))
                                
