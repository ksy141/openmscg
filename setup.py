import sys, os, re
sys.path.insert(0, '.')

from setuptools import setup, find_packages, Extension
import configparser
import numpy as np

package_name = 'mscg'

def read_version(version_file):
    with open(version_file, "r") as f:
        version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", f.read(), re.M)
    
    if version_match:
        return version_match.group(1)
    
    raise RuntimeError("Unable to find version string.")

def read_config(args):
    if not os.path.isfile('build.cfg'):
        return
    
    cfg = configparser.ConfigParser()
    cfg.readfp(open('build.cfg'))
    
    if not cfg.has_section('build_options'):
        return
    
    for key, val in cfg.items('build_options'):
        if key in args:
            if type(args[key]) == list:
                args[key] = val.strip().split()
            else:
                args[key] = val
        else:
            error = "Unknown key [" + key + "] in the configuration seciton [build_options]."
            raise RuntimeError(error)
    
    return

def parse_args(args):
    cli_args = sys.argv[1:]
    print("CLI Args: ", cli_args)
    
    for arg in cli_args:
        if arg.startswith('--with-'):
            w = arg[7:].split('=', 1)
            key = w[0].replace('-','_')
            
            if key in args:
                if type(args[key]) == list:
                    args[key] = w[1].split(' ') if len(w)>1 else []
                else:
                    args[key] = w[1] if len(w)>1 else ''
                    
                sys.argv.remove(arg)
    
def update_envs(args):
    if args['cc'] != '':
        os.environ['CC'] = args['cc']
    
    if setup_args['cxx'] != '':
        os.environ['CXX'] = args['cxx']

def read_requirements():
    with open("requirements.txt") as f:
        return f.read().strip().split("\n")

def build_defs(args):
    defs = []
    
    for arg in args:
        w = arg.split('=')
        defs.append((w[0], w[1]) if len(w)>1 else (w[0], None))
    
    return defs
    
    

if __name__ == '__main__':
    print("ENV=> ", os.environ)    
    
    setup_args = {
        'cc'        : '',
        'cxx'       : '',
        'compile'   : ['-O2', '-Wno-sign-compare'],
        'link'      : [],
        'gsl_lib'   : ['-lgsl', '-lgslcblas'],
        'lapack_def': ['USE_MKL'],
        'lapack_lib': ['-lmkl_gf_lp64', '-lmkl_intel_thread', '-lmkl_core'],
    }
    
    read_config(setup_args)
    parse_args(setup_args)
    update_envs(setup_args)
        
    for k in setup_args:
        print("Build option: " + k + "=" + str(setup_args[k]))
        
    core_prefix = package_name + '.core.cxx_'
    core_root   = package_name + '/core/'
    api_path    = core_root + 'api/'
    src_path    = core_root + 'src/'
    inc_path    = [src_path, np.get_include()]
    
    for prefix in ['CONDA_PREFIX', 'PREFIX']:
        if prefix in os.environ:
            inc_path.append(os.environ[prefix] + '/include')
    
    def src_files(api, files):
        return [api_path + api + '.cpp'] + [src_path + file + '.cpp' for file in files]
        
    def model_extention(name, extra_src = [], extra_link = []):
        return Extension(core_prefix + 'model_' + name,
            include_dirs = inc_path,
            sources = src_files('py_model_' + name, ['model', 'model_' + name] + extra_src),
            extra_compile_args = setup_args['compile'],
            extra_link_args = extra_link + setup_args['link'],
        )
    
    def get_models(model_path):
        models = []
        
        for f in [f for f in os.listdir(model_path) if f.endswith('.cpp') and f.startswith('model_')]:
            f = f.split('.')[0]
            src = src_files('py_' + f, [f, 'model'] + (['bspline'] if 'bspline' in f else []))
            print(f, "<=", src)
            
            models.append(Extension(
                core_prefix + f,
                include_dirs = inc_path,
                sources = src,
                extra_compile_args = setup_args['compile'],
                extra_link_args = setup_args['gsl_lib'] if 'bspline' in f else [] + setup_args['link'],
            ))
        
        return models
        
    
    extensions = [
        Extension(core_prefix + 'topol',
            include_dirs = inc_path,
            sources = src_files('py_topol', ['topology']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        ),

        Extension(core_prefix + 'traj',
            include_dirs = inc_path + [src_path + 'gmx', src_path + 'dcd'],
            sources = src_files('py_traj', ['traj','traj_lammps',
                                            'traj_trr','gmx/xdrfile','gmx/xdrfile_trr', 
                                            'traj_dcd','dcd/dcdplugin'
                                           ]),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        ),

        Extension(core_prefix + 'pairlist',
            include_dirs = inc_path,
            sources = src_files('py_pairlist', ['pair_list']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        ),
        
        Extension(core_prefix + 'bondlist',
            include_dirs = inc_path,
            sources = src_files('py_bondlist', ['bond_list']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        ),
        
        Extension(core_prefix + 'bspline',
            include_dirs = inc_path,
            sources = src_files('py_bspline', ['bspline']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['gsl_lib'] + setup_args['link'],
        ),

        Extension(core_prefix + 'tables',
            include_dirs = inc_path,
            sources = src_files('py_tables', ['tables']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        ),
        
        Extension(core_prefix + 'force',
            include_dirs = inc_path,
            sources = src_files('py_force', ['force']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        ),

        Extension(core_prefix + 'nb3b',
            include_dirs = inc_path,
            sources = src_files('py_nb3b', ['nb3b']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        )
    ] + get_models(src_path)
    
    entry_points = {"console_scripts": [
        f.split('.')[0] + "=" + package_name + ".cli." + f.split('.')[0] + ":main" \
        for f in os.listdir(package_name + "/cli") if f.startswith("cg") and f.endswith(".py")
    ]}
    
    setup(
        version = read_version(package_name + '/__init__.py'),
        python_requires = '>=3.6',
        install_requires = read_requirements(),
        packages = find_packages(),
        ext_modules = extensions,
        entry_points = entry_points
    )
