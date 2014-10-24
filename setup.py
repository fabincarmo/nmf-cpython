from distutils.core import setup, Extension

module1 = Extension('nmf', sources=['nmf.c'],
                    include_dirs=['/usr/include'],
                    libraries=['cblas','atlas']
                    )
setup(name = 'nmf',
        version='1.0',
        description='This is my package',
        ext_modules = [module1])

