from distutils.core import setup, Extension

module1 = Extension('mmeval',
                    sources = ['mmeval.c'],
                    extra_compile_args = ['-std=c99', '-O2'])

setup (name = 'MMEval',
       version = '1.0',
       description = 'Evaluates molecular mechanics energies',
       ext_modules = [module1])
