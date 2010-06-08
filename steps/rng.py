from . import steps_swig
import _steps_swig

### Now defunct mesh saving/loading tool ###
# from steps_swig import loadASCII, saveASCII


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

class RNG(steps_swig.RNG):
    """
    
    """
    def __init__(self, *args): 
        """
        Construction::
        
            g = steps.geom.Geom()
            
        Create a geometry container object.
            
        Arguments: 
            None
        """
        this = _steps_swig.new_RNG(*args)
        try: self.this.append(this)
        except: self.this = this
        # set Geom object to do all the cleaning-up
        self.thisown = True

def create(*args):
  """
    Creates and returns a reference to a steps.rng.RNG random number generator object, 
    which is specified by type and pre-allocates a buffer list with size of buffer_size.

    Syntax::
        
        create(type, buffer_size)

    Arguments:
        * string type
        * uint buffer_size

    Return:
        steps.rng.RNG

    """
  return steps_swig.create(*args)