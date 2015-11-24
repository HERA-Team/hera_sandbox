import glob

__version__ = '0.0.1'

setup_args = {
    'name': 'capo',
    'author': 'Danny Jacobs and Aaron Parsons',
    'author_email': 'jacobsda at sas.upenn.edu',
    'license': 'GPL',
    'package_dir' : {'capo':'src'},
    'packages' : ['capo'],
    'scripts': glob.glob('scripts/*'),
    'version': __version__,
    'package_data' : {'capo':['data/*.txt','data/psa32_apj_April15_zs_1/*.npz']}
}

if __name__ == '__main__':
    from distutils.core import setup
    apply(setup, (), setup_args)
