######################################################################
#
#   artemide for python3
#           A.Vladimirov (22.06.2020)
#
######################################################################

import artemide
import numpy

def initialize(fileName):
    """
    Initialization of artemide

    Parameters
    ----------
    fileName : string
        The path to the constants-file

    Returns
    -------
    None.

    """
    
    import os.path
    if os.path.exists(fileName):
    
        artemide.harpy.initialize(fileName)
        if artemide.harpy.started:
            pass
        else:
            print("Welcome to harpy -- the python interface for artemide")
    else:
        raise FileNotFoundError('consts-file '+fileName+' NOT FOUND')

def ShowStatistics():
    """
    Print the statistics

    Returns
    -------
    None.

    """
    artemide.harpy.showstatistics()

def setNPparameters(l):
    """
    Set NP parameters 

    Parameters
    ----------
    l : float array
        list of NP parameters

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_main(numpy.asfortranarray(l))
    else:
        raise TypeError()
        
def setNPparameters_TMDR(l):
    """
    Setting NP parameters for the model of TMDR

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_tmdr(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_tmdr(l)
    else:
        raise TypeError()
        

def setNPparameters_uTMDPDF(l):
    """
    Setting NP parameters for the model of uTMDPDF

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_utmdpdf(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_utmdpdf(l)
    else:
        raise TypeError()
        
        
def setNPparameters_uTMDFF(l):
    """
    Setting NP parameters for the model of uTMDFF

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_utmdff(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_utmdff(l)
    else:
        raise TypeError()
        
def setNPparameters_lpTMDPDF(l):
    """
    Setting NP parameters for the model of lpTMDPDF

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_lptmdpdf(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_lptmdpdf(l)
    else:
        raise TypeError()

def setNPparameters_SiversTMDPDF(l):
    """
    Setting NP parameters for the model of SiversTMDPDF

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_siverstmdpdf(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_siverstmdpdf(l)
    else:
        raise TypeError()
    
def setNPparameters_wgtTMDPDF(l):
    """
    Setting NP parameters for the model of wgtTMDPDF

    Parameters
    ----------
    l : float array or integer
        list of NP parameters, or the number of replica (if supported by model)

    Returns
    -------
    None.

    """
    if isinstance(l,list) or isinstance(l,numpy.ndarray):
        artemide.harpy.setlambda_wgttmdpdf(numpy.asfortranarray(l))
    elif isinstance(l, int):
        artemide.harpy.setreplica_wgttmdpdf(l)
    else:
        raise TypeError()


def varyScales(c1,c2,c3,c4):
    """
    Set new scale variation parameters

    Parameters
    ----------
    c1 : float
        Scale variation parameter c1
    c2 : float
        Scale variation parameter c2
    c3 : float
        Scale variation parameter c3
    c4 : float
        Scale variation parameter c4

    Returns
    -------
    None.

    """
    if not isinstance(c1,float):
        raise TypeError("c1 is not float")
    if not isinstance(c2,float):
        raise TypeError("c2 is not float")
    if not isinstance(c3,float):
        raise TypeError("c3 is not float")
    if not isinstance(c4,float):
        raise TypeError("c4 is not float")
        
    artemide.harpy.setscalevariation(c1,c2,c3,c4)

def _IsKinematicProper(s,qT,Q,y):
    """ Checks the point for the proper kinematics
    Especially for the correct domain of X.
    """
    gridX=0.00001
    if qT[0]>qT[1]:
        print('Wrong order of qT')
        return False
    if Q[0]>Q[1]:
        print('Wrong order of Q')
        return False
    if y[0]>y[1]:
        print('Wrong order of y')
        return False
    if qT[1]>Q[0]:
        print('qT (',qT[1],') > Q(', Q[0],')')
        return False
    if Q[1]>s:
        print('Q (',Q[1],') > s(', s,')')
        return False
    
    x1x2=(Q[0]**2+qT[0]**2)/s
    ymax=-numpy.log(numpy.sqrt(x1x2))
    ymin=-ymax
    
    if y[1]<ymin or y[0]>ymax:
        print('y ',y, 'is outside physical region ',[ymin,ymax])
        return False
    
    if y[1]>ymax:
        if x1x2<gridX:
            print('x outside of the grid')
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(-y[1])<gridX:
            print('x outside of the grid')
            return False
    
    if y[0]<ymin:
        if x1x2<gridX:
            print('x outside of the grid')
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(y[0])<gridX:
            print('x outside of the grid')
            return False
     
    x1x2=(Q[0]**2+qT[0]**2)/s
    ymax=-numpy.log(numpy.sqrt(x1x2))
    ymin=-ymax
    
    if y[1]<ymin or y[0]>ymax:
        print('y ',y, 'is outside physical region ',[ymin,ymax])
        return False
    
    if y[1]>ymax:
        if x1x2<gridX:
            print('x outside of the grid')
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(-y[1])<gridX:
            print('x outside of the grid')
            return False
    
    if y[0]<ymin:
        if x1x2<gridX:
            print('x outside of the grid')
            return False
    else:
        if numpy.sqrt(x1x2)*numpy.exp(y[0])<gridX:
            print('x outside of the grid')
            return False
    
    return True

def setPDFreplica(n,h=1):
    """
    Changes the replica for PDF input.
    
    This is a temporary function will be changed in future versions

    Parameters
    ----------
    n : Integer
        Number of PDF replica

    Returns
    -------
    None.

    """
    if not isinstance(n,int):
        raise TypeError()
    artemide.harpy.setpdfreplica(n,h)
    
def sethPDFreplica(n,h=1):
    """
    Changes the replica for PDF input.
    
    This is a temporary function will be changed in future versions

    Parameters
    ----------
    n : Integer
        Number of PDF replica

    Returns
    -------
    None.

    """
    if not isinstance(n,int):
        raise TypeError()
    artemide.harpy.setwgtpdfreplica(n,h)

def get_DNP(b,mu,f=1):
    
    if not isinstance(b,float):
        raise ValueError("parameter b must be float")
    elif (b<0.):
        b_internal=-b
    else:
        b_internal=b
    if not isinstance(mu,float):
        raise ValueError("parameter mu must be float")
    elif (mu<1.):
        raise ValueError("parameter mu must be > 1 GeV")
    if not isinstance(f, int):
        raise ValueError("parameter h must be integer")
        
    if(f==0):
        return artemide.harpy.getdnp(b_internal,mu,0)
    else:
        return artemide.harpy.getdnp(b_internal,mu,1)

def get_R(b,mu,zeta,f=1):

    if not isinstance(b,float):
        raise ValueError("parameter b must be float")
    elif (b<0.):
        b_internal=-b
    else:
        b_internal=b
    if not isinstance(mu,float):
        raise ValueError("parameter mu must be float")
    elif (mu<1.):
        raise ValueError("parameter mu must be > 1 GeV")

    if not isinstance(zeta,float):
        raise ValueError("parameter zeta must be float")
    elif (zeta<1.):
        raise ValueError("parameter zeta must be > 1 GeV")

    if not isinstance(f, int):
        raise ValueError("parameter h must be integer")


    if(f==0):
        return artemide.harpy.getr(b_internal,mu,zeta,0)
    else:
        return artemide.harpy.getr(b_internal,mu,zeta,1)

def get_uTMDPDF(x,b,h,mu=-1.,zeta=-1.,includeGluon=False):
    """
    Return the string of unpolarized TMDPDF 
    (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)
    If both mu and zeta are not specified (or negative) return optimal
    if mu is positive returns evolved to mu,mu^2
    
     

    Parameters
    ----------
    x : float in [0,1]
        Bjorken x
    b : float >0
        spatial parameter
    h : integer 1,2,3,...
        Hadron number
    mu : float, optional
        Scale of TMD mu [GeV]. The default is -1.
    zeta : float, optional
        Scale of TMD zeta [GeV]. The default is -1.
    includeGluon : bool, optional
        Include gluons or not. The default is False.

    Returns
    -------
    [list of float]
       (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)

    """
    
    if not isinstance(x, float):
        raise ValueError("parameter x must be float")
    elif (x<0.) or (x>1.):
        raise ValueError("parameter x must be in [0,1]")
    if not isinstance(b, float):
        raise ValueError("parameter b must be float")
    elif (b<0.):
        b_internal=-b
    else:
        b_internal=b
    if not isinstance(h, int):
        raise ValueError("parameter h must be integer")
    elif h<1:
        raise ValueError("parameter h expected to be positive integer")
        
    if not isinstance(mu, float):
        raise ValueError("parameter mu must be float")
    if not isinstance(zeta, float):
        raise ValueError("parameter zeta must be float")
    
    if mu<0.:
        if includeGluon:
            return artemide.harpy.utmdpdf_50_optimal(x,b_internal,h)
        else:
            return artemide.harpy.utmdpdf_5_optimal(x,b_internal,h)
    elif zeta<0:
        if includeGluon:
            return artemide.harpy.utmdpdf_50_evolved(x,b_internal,mu,mu**2,h)
        else:
            return artemide.harpy.utmdpdf_5_evolved(x,b_internal,mu,mu**2,h)
    else:
        if includeGluon:
            return artemide.harpy.utmdpdf_50_evolved(x,b_internal,mu,zeta,h)
        else:
            return artemide.harpy.utmdpdf_5_evolved(x,b_internal,mu,zeta,h)
        
def get_uTMDFF(x,b,h,mu=-1.,zeta=-1.,includeGluon=False):
    """
    Return the string of unpolarized TMDFF 
    (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)
    If both mu and zeta are not specified (or negative) return optimal
    if mu is positive returns evolved to mu,mu^2
    
     

    Parameters
    ----------
    x : float in [0,1]
        Bjorken x
    b : float >0
        spatial parameter
    h : integer 1,2,3,...
        Hadron number
    mu : float, optional
        Scale of TMD mu [GeV]. The default is -1.
    zeta : float, optional
        Scale of TMD zeta [GeV]. The default is -1.
    includeGluon : bool, optional
        Include gluons or not. The default is False.

    Returns
    -------
    [list of float]
       (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)

    """
    
    if not isinstance(x, float):
        raise ValueError("parameter x must be float")
    elif (x<0.) or (x>1.):
        raise ValueError("parameter x must be in [0,1]")
    if not isinstance(b, float):
        raise ValueError("parameter b must be float")
    elif (b<0.):
        b_internal=-b
    else:
        b_internal=b
    if not isinstance(h, int):
        raise ValueError("parameter x must be float")
    elif h<1:
        raise ValueError("parameter h expected to be positive integer")
        
    if not isinstance(mu, float):
        raise ValueError("parameter mu must be float")
    if not isinstance(zeta, float):
        raise ValueError("parameter zeta must be float")
    
    if mu<0.:
        if includeGluon:
            return artemide.harpy.utmdff_50_optimal(x,b_internal,h)
        else:
            return artemide.harpy.utmdff_5_optimal(x,b_internal,h)
    elif zeta<0:
        if includeGluon:
            return artemide.harpy.utmdff_50_evolved(x,b_internal,mu,mu**2,h)
        else:
            return artemide.harpy.utmdff_5_evolved(x,b_internal,mu,mu**2,h)
    else:
        if includeGluon:
            return artemide.harpy.utmdff_50_evolved(x,b_internal,mu,zeta,h)
        else:
            return artemide.harpy.utmdff_5_evolved(x,b_internal,mu,zeta,h)
        
def get_lpTMDPDF(x,b,h,mu=-1.,zeta=-1.,includeGluon=False):
    """
    Return the string of linearly polarized gluon TMDFF 
    (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)
    If both mu and zeta are not specified (or negative) return optimal
    if mu is positive returns evolved to mu,mu^2
    
     

    Parameters
    ----------
    x : float in [0,1]
        Bjorken x
    b : float >0
        spatial parameter
    h : integer 1,2,3,...
        Hadron number
    mu : float, optional
        Scale of TMD mu [GeV]. The default is -1.
    zeta : float, optional
        Scale of TMD zeta [GeV]. The default is -1.
    includeGluon : bool, optional
        IGNORED

    Returns
    -------
    [list of float]
       (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)

    """
    
    if not isinstance(x, float):
        raise ValueError("parameter x must be float")
    elif (x<0.) or (x>1.):
        raise ValueError("parameter x must be in [0,1]")
    if not isinstance(b, float):
        raise ValueError("parameter b must be float")
    elif (b<0.):
        b_internal=-b
    else:
        b_internal=b
    if not isinstance(h, int):
        raise ValueError("parameter x must be float")
    elif h<1:
        raise ValueError("parameter h expected to be positive integer")
        
    if not isinstance(mu, float):
        raise ValueError("parameter mu must be float")
    if not isinstance(zeta, float):
        raise ValueError("parameter zeta must be float")
    
    if mu<0.:
        return artemide.harpy.lptmdpdf_50_optimal(x,b_internal,h)
    elif zeta<0:
        return artemide.harpy.lptmdpdf_50_evolved(x,b_internal,mu,mu**2,h)        
    else:
        return artemide.harpy.lptmdpdf_50_evolved(x,b_internal,mu,zeta,h)
    
def get_SiversTMDPDF(x,b,h,mu=-1.,zeta=-1.,includeGluon=False):
    """
    Return the string of Sivers TMDPDF 
    (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)
    If both mu and zeta are not specified (or negative) return optimal
    if only mu is positive returns evolved to mu,mu^2
    if mu and zeta are positive returns evolved to mu,zeta
    
     

    Parameters
    ----------
    x : float in [0,1]
        Bjorken x
    b : float >0
        spatial parameter
    h : integer 1,2,3,...
        Hadron number
    mu : float, optional
        Scale of TMD mu [GeV]. The default is -1.
    zeta : float, optional
        Scale of TMD zeta [GeV]. The default is -1.
    includeGluon : bool, optional
        Include gluons or not. The default is False.

    Returns
    -------
    [list of float]
       (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)

    """
    
    if not isinstance(x, float):
        raise ValueError("parameter x must be float")
    elif (x<0.) or (x>1.):
        raise ValueError("parameter x must be in [0,1]")
    if not isinstance(b, float):
        raise ValueError("parameter b must be float")
    elif (b<0.):
        b_internal=-b
    else:
        b_internal=b
    if not isinstance(h, int):
        raise ValueError("parameter x must be float")
    elif h<1:
        raise ValueError("parameter h expected to be positive integer")
        
    if not isinstance(mu, float):
        raise ValueError("parameter mu must be float")
    if not isinstance(zeta, float):
        raise ValueError("parameter zeta must be float")
    
    if mu<0.:
        if includeGluon:
            return artemide.harpy.siverstmdpdf_50_optimal(x,b_internal,h)
        else:
            return artemide.harpy.siverstmdpdf_5_optimal(x,b_internal,h)
    elif zeta<0:
        if includeGluon:
            return artemide.harpy.siverstmdpdf_50_evolved(x,b_internal,mu,mu**2,h)
        else:
            return artemide.harpy.siverstmdpdf_5_evolved(x,b_internal,mu,mu**2,h)
    else:
        if includeGluon:
            return artemide.harpy.siverstmdpdf_50_evolved(x,b_internal,mu,zeta,h)
        else:
            return artemide.harpy.siverstmdpdf_5_evolved(x,b_internal,mu,zeta,h)
        
def get_wgtTMDPDF(x,b,h,mu=-1.,zeta=-1.,includeGluon=False):
    """
    Return the string of wgt TMDPDF 
    (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)
    If both mu and zeta are not specified (or negative) return optimal
    if only mu is positive returns evolved to mu,mu^2
    if mu and zeta are positive returns evolved to mu,zeta
    
     

    Parameters
    ----------
    x : float in [0,1]
        Bjorken x
    b : float >0
        spatial parameter
    h : integer 1,2,3,...
        Hadron number
    mu : float, optional
        Scale of TMD mu [GeV]. The default is -1.
    zeta : float, optional
        Scale of TMD zeta [GeV]. The default is -1.
    includeGluon : bool, optional
        Include gluons or not. The default is False.

    Returns
    -------
    [list of float]
       (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)

    """
    
    if not isinstance(x, float):
        raise ValueError("parameter x must be float")
    elif (x<0.) or (x>1.):
        raise ValueError("parameter x must be in [0,1]")
    if not isinstance(b, float):
        raise ValueError("parameter b must be float")
    elif (b<0.):
        b_internal=-b
    else:
        b_internal=b
    if not isinstance(h, int):
        raise ValueError("parameter x must be float")
    elif h<1:
        raise ValueError("parameter h expected to be positive integer")
        
    if not isinstance(mu, float):
        raise ValueError("parameter mu must be float")
    if not isinstance(zeta, float):
        raise ValueError("parameter zeta must be float")
    
    if mu<0.:
        if includeGluon:
            return artemide.harpy.wgttmdpdf_50_optimal(x,b_internal,h)
        else:
            return artemide.harpy.wgttmdpdf_5_optimal(x,b_internal,h)
    elif zeta<0:
        if includeGluon:
            return artemide.harpy.wgttmdpdf_50_evolved(x,b_internal,mu,mu**2,h)
        else:
            return artemide.harpy.wgttmdpdf_5_evolved(x,b_internal,mu,mu**2,h)
    else:
        if includeGluon:
            return artemide.harpy.wgttmdpdf_50_evolved(x,b_internal,mu,zeta,h)
        else:
            return artemide.harpy.wgttmdpdf_5_evolved(x,b_internal,mu,zeta,h)

###############################################################################

def get_uTMDPDF_kT(x,kT,h,mu=-1.,zeta=-1.,includeGluon=False):
    """
    Return the string of unpolarized TMDPDF in kT-space
    (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)
    If both mu and zeta are not specified (or negative) return optimal
    if mu is positive returns evolved to mu,mu^2
    
     

    Parameters
    ----------
    x : float in [0,1]
        Bjorken x
    kT : float >0
        spatial parameter
    h : integer 1,2,3,...
        Hadron number
    mu : float, optional
        Scale of TMD mu [GeV]. The default is -1.
    zeta : float, optional
        Scale of TMD zeta [GeV]. The default is -1.
    includeGluon : bool, optional
        Include gluons or not. The default is False.

    Returns
    -------
    [list of float]
       (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)

    """
    
    if not isinstance(x, float):
        raise ValueError("parameter x must be float")
    elif (x<0.) or (x>1.):
        raise ValueError("parameter x must be in [0,1]")
    if not isinstance(kT, float):
        raise ValueError("parameter kT must be float")
    elif (kT<0.):
        kT_internal=-kT
    else:
        kT_internal=kT
    if not isinstance(h, int):
        raise ValueError("parameter h must be integer")
    elif h<1:
        raise ValueError("parameter h expected to be positive integer")
        
    if not isinstance(mu, float):
        raise ValueError("parameter mu must be float")
    if not isinstance(zeta, float):
        raise ValueError("parameter zeta must be float")
    
    if mu<0.:
        if includeGluon:
            return artemide.harpy.utmdpdf_kt_50_optimal(x,kT_internal,h)
        else:
            return artemide.harpy.utmdpdf_kt_5_optimal(x,kT_internal,h)
    elif zeta<0:
        if includeGluon:
            return artemide.harpy.utmdpdf_kt_50_evolved(x,kT_internal,mu,mu**2,h)
        else:
            return artemide.harpy.utmdpdf_kt_5_evolved(x,kT_internal,mu,mu**2,h)
    else:
        if includeGluon:
            return artemide.harpy.utmdpdf_kt_50_evolved(x,kT_internal,mu,zeta,h)
        else:
            return artemide.harpy.utmdpdf_kt_5_evolved(x,kT_internal,mu,zeta,h)
        
def get_uTMDFF_kT(x,kT,h,mu=-1.,zeta=-1.,includeGluon=False):
    """
    Return the string of unpolarized TMDFF  in kT-space
    (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)
    If both mu and zeta are not specified (or negative) return optimal
    if mu is positive returns evolved to mu,mu^2
    
     

    Parameters
    ----------
    x : float in [0,1]
        Bjorken x
    kT : float >0
        spatial parameter
    h : integer 1,2,3,...
        Hadron number
    mu : float, optional
        Scale of TMD mu [GeV]. The default is -1.
    zeta : float, optional
        Scale of TMD zeta [GeV]. The default is -1.
    includeGluon : bool, optional
        Include gluons or not. The default is False.

    Returns
    -------
    [list of float]
       (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)

    """
    
    if not isinstance(x, float):
        raise ValueError("parameter x must be float")
    elif (x<0.) or (x>1.):
        raise ValueError("parameter x must be in [0,1]")
    if not isinstance(kT, float):
        raise ValueError("parameter kT must be float")
    elif (kT<0.):
        kT_internal=-kT
    else:
        kT_internal=kT
    if not isinstance(h, int):
        raise ValueError("parameter x must be float")
    elif h<1:
        raise ValueError("parameter h expected to be positive integer")
        
    if not isinstance(mu, float):
        raise ValueError("parameter mu must be float")
    if not isinstance(zeta, float):
        raise ValueError("parameter zeta must be float")
    
    if mu<0.:
        if includeGluon:
            return artemide.harpy.utmdff_kt_50_optimal(x,kT_internal,h)
        else:
            return artemide.harpy.utmdff_kt_5_optimal(x,kT_internal,h)
    elif zeta<0:
        if includeGluon:
            return artemide.harpy.utmdff_kt_50_evolved(x,kT_internal,mu,mu**2,h)
        else:
            return artemide.harpy.utmdff_kt_5_evolved(x,kT_internal,mu,mu**2,h)
    else:
        if includeGluon:
            return artemide.harpy.utmdff_kt_50_evolved(x,kT_internal,mu,zeta,h)
        else:
            return artemide.harpy.utmdff_kt_5_evolved(x,kT_internal,mu,zeta,h)
        
def get_lpTMDPDF_kT(x,kT,h,mu=-1.,zeta=-1.,includeGluon=False):
    """
    Return the string of linearly polarized gluon TMDFF  in kT-space
    (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)
    If both mu and zeta are not specified (or negative) return optimal
    if mu is positive returns evolved to mu,mu^2
    
     

    Parameters
    ----------
    x : float in [0,1]
        Bjorken x
    kT : float >0
        spatial parameter
    h : integer 1,2,3,...
        Hadron number
    mu : float, optional
        Scale of TMD mu [GeV]. The default is -1.
    zeta : float, optional
        Scale of TMD zeta [GeV]. The default is -1.
    includeGluon : bool, optional
        IGNORED

    Returns
    -------
    [list of float]
       (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)

    """
    
    if not isinstance(x, float):
        raise ValueError("parameter x must be float")
    elif (x<0.) or (x>1.):
        raise ValueError("parameter x must be in [0,1]")
    if not isinstance(kT, float):
        raise ValueError("parameter kT must be float")
    elif (kT<0.):
        kT_internal=-kT
    else:
        kT_internal=kT
    if not isinstance(h, int):
        raise ValueError("parameter x must be float")
    elif h<1:
        raise ValueError("parameter h expected to be positive integer")
        
    if not isinstance(mu, float):
        raise ValueError("parameter mu must be float")
    if not isinstance(zeta, float):
        raise ValueError("parameter zeta must be float")
    
    if mu<0.:
        return artemide.harpy.lptmdpdf_kt_50_optimal(x,kT_internal,h)
    elif zeta<0:
        return artemide.harpy.lptmdpdf_kt_50_evolved(x,kT_internal,mu,mu**2,h)        
    else:
        return artemide.harpy.lptmdpdf_kt_50_evolved(x,kT_internal,mu,zeta,h)
    
def get_SiversTMDPDF_kT(x,kT,h,mu=-1.,zeta=-1.,includeGluon=False):
    """
    Return the string of Sivers TMDPDF  in kT-space
    (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)
    If both mu and zeta are not specified (or negative) return optimal
    if mu is positive returns evolved to mu,mu^2
    
     

    Parameters
    ----------
    x : float in [0,1]
        Bjorken x
    kT : float >0
        spatial parameter
    h : integer 1,2,3,...
        Hadron number
    mu : float, optional
        Scale of TMD mu [GeV]. The default is -1.
    zeta : float, optional
        Scale of TMD zeta [GeV]. The default is -1.
    includeGluon : bool, optional
        Include gluons or not. The default is False.

    Returns
    -------
    [list of float]
       (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)

    """
    
    if not isinstance(x, float):
        raise ValueError("parameter x must be float")
    elif (x<0.) or (x>1.):
        raise ValueError("parameter x must be in [0,1]")
    if not isinstance(kT, float):
        raise ValueError("parameter kT must be float")
    elif (kT<0.):
        kT_internal=-kT
    else:
        kT_internal=kT
    if not isinstance(h, int):
        raise ValueError("parameter x must be float")
    elif h<1:
        raise ValueError("parameter h expected to be positive integer")
        
    if not isinstance(mu, float):
        raise ValueError("parameter mu must be float")
    if not isinstance(zeta, float):
        raise ValueError("parameter zeta must be float")
    
    if mu<0.:
        if includeGluon:
            return artemide.harpy.siverstmdpdf_kt_50_optimal(x,kT_internal,h)
        else:
            return artemide.harpy.siverstmdpdf_kt_5_optimal(x,kT_internal,h)
    elif zeta<0:
        if includeGluon:
            return artemide.harpy.siverstmdpdf_kt_50_evolved(x,kT_internal,mu,mu**2,h)
        else:
            return artemide.harpy.siverstmdpdf_kt_5_evolved(x,kT_internal,mu,mu**2,h)
    else:
        if includeGluon:
            return artemide.harpy.siverstmdpdf_kt_50_evolved(x,kT_internal,mu,zeta,h)
        else:
            return artemide.harpy.siverstmdpdf_kt_5_evolved(x,kT_internal,mu,zeta,h)
        
def get_wgtTMDPDF_kT(x,kT,h,mu=-1.,zeta=-1.,includeGluon=False):
    """
    Return the string of Worm-gear T TMDPDF  in kT-space
    (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)
    If both mu and zeta are not specified (or negative) return optimal
    if mu is positive returns evolved to mu,mu^2
    
     

    Parameters
    ----------
    x : float in [0,1]
        Bjorken x
    kT : float >0
        spatial parameter
    h : integer 1,2,3,...
        Hadron number
    mu : float, optional
        Scale of TMD mu [GeV]. The default is -1.
    zeta : float, optional
        Scale of TMD zeta [GeV]. The default is -1.
    includeGluon : bool, optional
        Include gluons or not. The default is False.

    Returns
    -------
    [list of float]
       (bbar,cbar,sbar,ubar,dbar,gluon,d,u,s,c,b)

    """
    
    if not isinstance(x, float):
        raise ValueError("parameter x must be float")
    elif (x<0.) or (x>1.):
        raise ValueError("parameter x must be in [0,1]")
    if not isinstance(kT, float):
        raise ValueError("parameter kT must be float")
    elif (kT<0.):
        kT_internal=-kT
    else:
        kT_internal=kT
    if not isinstance(h, int):
        raise ValueError("parameter x must be float")
    elif h<1:
        raise ValueError("parameter h expected to be positive integer")
        
    if not isinstance(mu, float):
        raise ValueError("parameter mu must be float")
    if not isinstance(zeta, float):
        raise ValueError("parameter zeta must be float")
    
    if mu<0.:
        if includeGluon:
            return artemide.harpy.wgttmdpdf_kt_50_optimal(x,kT_internal,h)
        else:
            return artemide.harpy.wgttmdpdf_kt_5_optimal(x,kT_internal,h)
    elif zeta<0:
        if includeGluon:
            return artemide.harpy.wgttmdpdf_kt_50_evolved(x,kT_internal,mu,mu**2,h)
        else:
            return artemide.harpy.wgttmdpdf_kt_5_evolved(x,kT_internal,mu,mu**2,h)
    else:
        if includeGluon:
            return artemide.harpy.wgttmdpdf_kt_50_evolved(x,kT_internal,mu,zeta,h)
        else:
            return artemide.harpy.wgttmdpdf_kt_5_evolved(x,kT_internal,mu,zeta,h)


###############################################################################
class DY:
        """Static class for evaluation of DY cross-section
        """

        @staticmethod
        def xSec(process,s,qT,Q,y,includeCuts,CutParameters=None,Num=None):
                """Cross-section for DY integrated over bin
                      
                Arguments: (process,s,qT,Q,y,includeCuts,CutParameters=None,Num=4)
                process         = (int, int, int) (see definition in artemide manual)
                s               = Mandelshtam variable s
                qT              = (qT-Min,qT-Max) boundaries of bin in qT
                Q               = (Q-Min,Q-Max) boundaries of bin in Q
                y               = (y-Min,y-Max) boundaries of bin in y
                includeCuts     = True/False to include leptonic cuts
                CutParameters   = (real,real,real,real) must be if includeCuts=True (see definition in artemide manual)
                Num             = even integer, number of section of qt-integration (defaul=4)
                """

                
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
                                
