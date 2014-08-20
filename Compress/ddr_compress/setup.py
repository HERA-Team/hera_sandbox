import glob
__version__ = '0.0.1'

setup_args = {
    'name': 'ddr_compress',
    'author': 'PAPER team',
    'author_email': 'jacobsda at sas.upenn.edu',
    'license': 'GPL',
    'package_dir' : {'ddr_compress':'src'},
    'packages' : ['ddr_compress'],
    'scripts': glob.glob('scripts/*'),
    'version': __version__,
}

if __name__ == '__main__':
    from distutils.core import setup
    apply(setup, (), setup_args)
