from argschema import ArgSchema, ArgSchemaParser 
from argschema.schemas import DefaultSchema
from argschema.fields import Nested, InputDir, String, Float, Dict, Int, Bool, NumpyArray
from ...common.schemas import EphysParams, Directories, CommonFiles



class KS4HelperParameters(DefaultSchema):

    do_CAR = Bool(required=False, default=True, help='set to True to perform common average referencing (median subtraction)')
    Th_universal = Float(required=False, default=9, help='threshold for creating templates')
    Th_learned = Float(required=False, default=8, help='threshold for creating templates')
    Th_single_chan = Float(required=False, default=8, help='threshold crossings for pre-clustering (in PCA projection space)')
    duplicate_spike_bins = Int(required=False, default=15, help='Number of bins for which subsequent spikes from the same cluster are assumed to be artifacts. A value of 0 disables this step.')
    nblocks = Int(required=False, default=5, help='number of blocks used to segment the probe when tracking drift, 0 == do not track, 1 == rigid, > 1 == non-rigid')
    sig_interp = Float(required=True, default=20.0, help='sigma for the Gaussian interpolation in drift correction um)')
    whitening_range = Int(required=False, default=32, help='number of channels to use for whitening each channel')
    min_template_size = Float(required=False, default=10, help='Width in um of Gaussian envelope for template weight')
    template_sizes = Int(required=False, default=5, help='number of template sizes, multiples of min_template size')
    template_from_data=Bool(required=False, default=True, help='set to True to extract templates from data')
    tmin = Float(required=False, default=0, help='time in sec to start processing')
    tmax = Float(required=False, default=-1, help='time in sec to end processing; if < 0, set to inf ')
    nearest_chans = Int(required=False, default=10, help='Number of nearest channels to consider when finding local maxima during spike detection.')
    nearest_templates = Int(required=False, default=100, help='Number of nearest spike template locations to consider when finding local maxima during spike detection.')
    ccg_threshold = Float(required=False, default=0.25, help='Fraction of refractory period violations that are allowed in the CCG compared to baseline; used to perform splits and merges. ')
    acg_threshold = Float(required=False, default=0.20, help='Fraction of refractory period violations that are allowed in the ACG compared to baseline; used to assign "good" units. ')
    ks_make_copy = Bool(required=False, default=False, help='If true, make a copy of the original KS output')
    
    doFilter = Int(required=False, default=0, help='filter if = 1, skip bp filtering if = 0')
    templateSeed = Int(required=False, default=0, help='filter if = 1, seed to pick which batches are used for finding universal templates')

class InputParameters(ArgSchema):
    ks4_helper_params = Nested(KS4HelperParameters)
    directories = Nested(Directories)
    ephys_params = Nested(EphysParams)
    common_files = Nested(CommonFiles)
    

class OutputSchema(DefaultSchema): 
    input_parameters = Nested(InputParameters, 
                              description=("Input parameters the module " 
                                           "was run with"), 
                              required=True) 
 
 
class OutputParameters(OutputSchema): 
    message = String()
    execution_time = Float()
