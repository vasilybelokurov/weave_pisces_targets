import os.path
import yaml
import astropy.table as atpy
import numpy as np

dirpref = os.path.dirname(__file__)

confdir = os.path.dirname(__file__) + '/../conf/'
with open(confdir + '/target.yaml', 'r') as fp:
    conf = yaml.load(fp, Loader=yaml.FullLoader)
with open(confdir + '/general.yaml', 'r') as fp:
    confGeneral = yaml.load(fp, Loader=yaml.FullLoader)


def get_aps_flag_dict():
    # read the aps table into dictionary
    # and turn the dictionary that maps our internal object classes
    # into integer array [0,1,0,3,]
    aps_tab_dict = dict(
        zip(*np.loadtxt(dirpref + '/../data/aps_flag.dat', dtype=str).T))
    for k in aps_tab_dict.keys():
        arr = np.array(list(aps_tab_dict[k]), dtype=int)
        assert (arr.max() == 1)
        assert (arr.min() == 0)
        aps_tab_dict[k] = arr
    # aps_tab_dict maps the APS object type to mask
    # we need to convert it to our internal type
    retD = {}
    for k in get_object_types():
        retD[k] = aps_tab_dict[conf['aps_flag_map'][k]]
    retD['__STAR__'] = aps_tab_dict['STAR']
    # since I have to add a general STAR APS_FLAG category
    # star APS flag, I'm creating a fake entry here
    return retD


def get_wsdb_host():
    home = os.environ['HOME']
    wsdb_file = home + '/wsdb_host'
    if os.path.exists(wsdb_file):
        with open(wsdb_file, 'r') as fp:
            wsdb = fp.read()
    else:
        try:
            wsdb = os.environ['WSDB_HOST']
        except KeyError:
            print('''you must specify the host name of the database either
with the wsdb_host file in your $HOME directory or throguh a WSDB_HOST
environmental variable
            ''')
            raise
    return wsdb


def get_object_types():
    return conf['object_bitmask'].keys()


def get_priority(name):
    """ return priority range for a given target class"""
    if name not in conf['priority']:
        raise Exception('Unknown object type')
    return conf['priority'][name]


def get_object_bitmask(name):
    """ return bitmask value for a given target class"""
    if name not in conf['object_bitmask']:
        raise Exception('Unknown object type')
    return 1 << conf['object_bitmask'][name]


def get_object_bitmask_dict():
    """ return priority range for a given target class"""
    D = {}
    for name in conf['object_bitmask'].keys():
        D[name] = 1 << conf['object_bitmask'][name]
    return D


def get_object_targprog_dict():
    """ return priority range for a given target class"""
    D = {}
    for name in conf['targprog'].keys():
        D[name] = conf['targprog'][name]
    return D


def get_sky_density():
    return confGeneral['sky_density']
