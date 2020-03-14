import sys, os, re
sys.path.insert(0, '.')

from setuptools import setup, find_packages, Extension
import configparser

def read_version(version_file):
    
    with open(version_file, "r") as f:
        version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", f.read(), re.M)
    
    if version_match:
        return version_match.group(1)
    
    raise RuntimeError("Unable to find version string.")

def read_config(args):
    
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

def read_requirements():
    with open("requirements.txt") as f:
        return f.read().strip().split("\n")

if __name__ == '__main__':
    
    setup_args = {
        'cc'        : '',
        'cxx'       : '',
        'compile'   : [],
        'link'      : [],
        'gsl_lib'   : [],
        'lapack_lib': [],
    }
    
    read_config(setup_args)
    
    for k in setup_args:
        print("Build option: " + k + "=" + str(setup_args[k]))
    
    if setup_args['cc'] != '':
        os.environ['CC'] = setup_args['cc']
    
    if setup_args['cxx'] != '':
        os.environ['CXX'] = setup_args['cxx']
        
    core_prefix = 'cg.core.cxx_'
    core_root   = 'cg/core/'
    api_path    = core_root + 'api/'
    src_path    = core_root + 'src/'

    def src_files(api, files):
        return [api_path + api + '.cpp'] + [src_path + file + '.cpp' for file in files]
    
    def table_extention(name):
        return Extension(core_prefix + 'table_' + name + '_bspline',
            include_dirs = [src_path],
            sources = src_files('py_table_' + name + '_bspline', ['table', 'bspline', 'table_' + name + '_bspline']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['gsl_lib'] + setup_args['link'],
        )
    
    extensions = [
        Extension(core_prefix + 'topol',
            include_dirs = [src_path],
            sources = src_files('py_topol', ['topology']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        ),

        Extension(core_prefix + 'traj',
            include_dirs = [src_path],
            sources = src_files('py_traj', ['traj','traj_lammps','traj_trr','xdrfile','xdrfile_trr']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        ),

        Extension(core_prefix + 'pairlist',
            include_dirs = [src_path],
            sources = src_files('py_pairlist', ['pair_list']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        ),
        
        Extension(core_prefix + 'bondlist',
            include_dirs = [src_path],
            sources = src_files('py_bondlist', ['bond_list']),
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['link'],
        ),
        
        Extension(core_prefix + 'matrix',
            include_dirs = [src_path],
            sources = src_files('py_matrix', ['matrix']),
            define_macros=[('USE_MKL', None)],
            extra_compile_args = setup_args['compile'],
            extra_link_args = setup_args['lapack_lib'] + setup_args['link'],
        ),
        
        table_extention('pair'),
        table_extention('bond'),
        table_extention('angle')
    ]
    
    entry_points = {"console_scripts": [
        f.split('.')[0] + "=cg.cli." + f.split('.')[0] + ":main" \
        for f in os.listdir("cg/cli") if f.startswith("cg") and f.endswith(".py")
    ]}
    
    setup(
        version = read_version('cg/__init__.py'),
        python_requires = '>=3.6',
        install_requires = read_requirements(),
        packages = find_packages(),
        ext_modules = extensions,
        entry_points = entry_points
    )
