# model.py

"""
This file is the user-interface file for all model objects in steps.
All objects are directly derived from the corresponding swig objects.
Model container object is owned by Python
All other objects are owned by c++ and Model container is responsible for 
all the cleaning-up of these objects (see cpp/steps/model/model.cpp class destructor).

"""

import model_swig
import _model_swig


###########################################
"""
class kcst(object):
    def __init__(self, fget, fset, _fget, _fset) :
        self.fget = fget
        self.fset = fset
        self._fget = _fget
        self._fset = _fset
    def __get__(self, instance, cls):
        k = self._fget(instance)
        self.fset(instance,k)
        return self.fget(instance)
    def __set__(self, instance, value):
        __swig_setmethods__(instance, value)
        self.fset(instance, value)
    def __swig_getmethods__ (self,instance,cls):
        return self._fget(instance)
    def __swig_setmethods__(self,instance,value):
        self._fset(instance,value)
"""
###########################################

class Model(model_swig.Model) :
    def __init__(self, *args): 
        """__init__(self) -> Model"""
        this = _model_swig.new_Model(*args)
        try: self.this.append(this)
        except: self.this = this
        # let the model object do the cleaning-up
        self.thisown = True

class Spec(model_swig.Spec) :
    def __init__(self, *args): 
        """__init__(self, string id, Model model) -> Spec"""
        this = _model_swig.new_Spec(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.Spec_setID
        self.__swig_getmethods__["id"] = _model_swig.Spec_getID
        self.__swig_getmethods__["model"] = _model_swig.Spec_getModel
    id = model_swig._swig_property(_model_swig.Spec_getID, _model_swig.Spec_setID)
    model = model_swig._swig_property(_model_swig.Spec_getModel)

class Surfsys(model_swig.Surfsys):
    def __init__(self, *args): 
        """__init__(self, string id, Model model) -> Surfsys"""
        this = _model_swig.new_Surfsys(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.Surfsys_setID
        self.__swig_getmethods__["id"] = _model_swig.Surfsys_getID
        self.__swig_getmethods__["model"] = _model_swig.Surfsys_getModel
    id = model_swig._swig_property(_model_swig.Surfsys_getID, _model_swig.Surfsys_setID)
    model = model_swig._swig_property(_model_swig.Surfsys_getModel)

class Volsys(model_swig.Volsys) :
    def __init__(self, *args): 
        """__init__(self, string id, Model model) -> Volsys"""
        this = _model_swig.new_Volsys(*args)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.Volsys_setID
        self.__swig_getmethods__["id"] = _model_swig.Volsys_getID
        self.__swig_getmethods__["model"] = _model_swig.Volsys_getModel
    id = model_swig._swig_property(_model_swig.Volsys_getID, _model_swig.Volsys_setID)
    model = model_swig._swig_property(_model_swig.Volsys_getModel)

class Diff(model_swig.Diff) :
    def __init__(self, *args, **kwargs): 
        """__init__(self, string id, Volsys volsys, Spec lig, double dcst=0.0) -> Diff"""
        this = _model_swig.new_Diff(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.Diff_setID
        self.__swig_getmethods__["id"] = _model_swig.Diff_getID
        self.__swig_getmethods__["model"] = _model_swig.Diff_getModel
        self.__swig_getmethods__["volsys"] = _model_swig.Diff_getVolsys
        self.__swig_setmethods__["lig"] = _model_swig.Diff_setLig
        self.__swig_getmethods__["lig"] = _model_swig.Diff_getLig
        self.__swig_setmethods__["dcst"] = _model_swig.Diff_setDcst
        self.__swig_getmethods__["dcst"] = _model_swig.Diff_getDcst
    id = model_swig._swig_property(_model_swig.Diff_getID, _model_swig.Diff_setID)
    model = model_swig._swig_property(_model_swig.Diff_getModel)
    volsys = model_swig._swig_property(_model_swig.Diff_getVolsys)
    dcst = model_swig._swig_property(_model_swig.Diff_getDcst, _model_swig.Diff_setDcst)
    lig = model_swig._swig_property(_model_swig.Diff_getLig, _model_swig.Diff_setLig)

class Reac(model_swig.Reac) :
    def __init__(self, *args, **kwargs): 
        
        #__init__(self, string id, Volsys volsys, vector_spc lhs=std::vector< steps::model::Spec * >(), 
        #    vector_spc rhs=std::vector< steps::model::Spec * >(), 
        #    double kcst=0.0) -> Reac
        
        this = _model_swig.new_Reac(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.Reac_setID
        self.__swig_getmethods__["id"] = _model_swig.Reac_getID
        self.__swig_getmethods__["model"] = _model_swig.Reac_getModel
        self.__swig_getmethods__["volsys"] = _model_swig.Reac_getVolsys
        self.__swig_setmethods__["kcst"] = _model_swig.Reac_setKcst
        self.__swig_getmethods__["kcst"] = _model_swig.Reac_getKcst
        self.__swig_setmethods__["lhs"] = _model_swig.Reac_setLHS
        self.__swig_getmethods__["lhs"] = _model_swig.Reac_getLHS
        self.__swig_setmethods__["rhs"] = _model_swig.Reac_setRHS
        self.__swig_getmethods__["rhs"] = _model_swig.Reac_getRHS
        self.__swig_getmethods__["order"] = _model_swig.Reac_getOrder
    id = model_swig._swig_property(_model_swig.Reac_getID, _model_swig.Reac_setID)
    model = model_swig._swig_property(_model_swig.Reac_getModel)
    volsys = model_swig._swig_property(_model_swig.Reac_getVolsys)
    kcst = model_swig._swig_property(_model_swig.Reac_getKcst, _model_swig.Reac_setKcst)
    lhs = model_swig._swig_property(_model_swig.Reac_getLHS, _model_swig.Reac_setLHS)
    rhs = model_swig._swig_property(_model_swig.Reac_getRHS, _model_swig.Reac_setRHS)
    order = model_swig._swig_property(_model_swig.Reac_getOrder)

class SReac(model_swig.SReac) :
    def __init__(self, *args, **kwargs): 
        """
        __init__(self, string id, Surfsys surfsys, vector_spc vlhs=std::vector< steps::model::Spec * >(), 
            vector_spc slhs=std::vector< steps::model::Spec * >(), 
            vector_spc irhs=std::vector< steps::model::Spec * >(), 
            vector_spc srhs=std::vector< steps::model::Spec * >(), 
            vector_spc orhs=std::vector< steps::model::Spec * >(), 
            double kcst=0.0) -> SReac
        """
        this = _model_swig.new_SReac(*args, **kwargs)
        try: self.this.append(this)
        except: self.this = this
        self.thisown = False
        self.__swig_setmethods__["id"] = _model_swig.SReac_setID
        self.__swig_getmethods__["id"] = _model_swig.SReac_getID
        self.__swig_getmethods__["model"] = _model_swig.SReac_getModel
        self.__swig_getmethods__["surfsys"] = _model_swig.SReac_getSurfsys
        self.__swig_setmethods__["kcst"] = _model_swig.SReac_setKcst
        self.__swig_getmethods__["kcst"] = _model_swig.SReac_getKcst
        self.__swig_getmethods__["order"] = _model_swig.SReac_getOrder
        self.__swig_setmethods__["outer"] = _model_swig.SReac_setOuter
        self.__swig_getmethods__["outer"] = _model_swig.SReac_getOuter
        self.__swig_setmethods__["vlhs"] = _model_swig.SReac_setVLHS
        self.__swig_getmethods__["vlhs"] = _model_swig.SReac_getVLHS
        self.__swig_setmethods__["slhs"] = _model_swig.SReac_setSLHS
        self.__swig_getmethods__["slhs"] = _model_swig.SReac_getSLHS
        self.__swig_setmethods__["irhs"] = _model_swig.SReac_setIRHS
        self.__swig_getmethods__["irhs"] = _model_swig.SReac_getIRHS
        self.__swig_setmethods__["srhs"] = _model_swig.SReac_setSRHS
        self.__swig_getmethods__["srhs"] = _model_swig.SReac_getSRHS
        self.__swig_setmethods__["orhs"] = _model_swig.SReac_setORHS
        self.__swig_getmethods__["orhs"] = _model_swig.SReac_getORHS
    id = model_swig._swig_property(_model_swig.SReac_getID, _model_swig.SReac_setID)
    model = model_swig._swig_property(_model_swig.SReac_getModel)
    surfsys = model_swig._swig_property(_model_swig.SReac_getSurfsys)
    kcst = model_swig._swig_property(_model_swig.SReac_getKcst, _model_swig.SReac_setKcst)
    outer = model_swig._swig_property(_model_swig.SReac_getOuter, _model_swig.SReac_setOuter)
    vlhs = model_swig._swig_property(_model_swig.SReac_getVLHS, _model_swig.SReac_setVLHS)
    slhs = model_swig._swig_property(_model_swig.SReac_getSLHS, _model_swig.SReac_setSLHS)
    irhs = model_swig._swig_property(_model_swig.SReac_getIRHS, _model_swig.SReac_setIRHS)
    srhs = model_swig._swig_property(_model_swig.SReac_getSRHS, _model_swig.SReac_setSRHS)
    orhs = model_swig._swig_property(_model_swig.SReac_getORHS, _model_swig.SReac_setORHS)
    order = model_swig._swig_property(_model_swig.SReac_getOrder)



