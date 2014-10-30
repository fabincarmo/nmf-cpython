from distutils.core import setup, Extension

module1 = Extension('nmf', sources=['pynmf.c'],
                    include_dirs=['/usr/include'],
                    libraries=['cblas']
                    )
setup(name = 'nmf',
        version='1.0',
        description='This is my package',
        ext_modules = [module1])

