from setuptools import setup

setup(name='zeus',
      version='0.1',
      description='Zeus - Different exons splicing comparison',
      url='http://github.com/ptdtan/BENTO-Seq/zeus',
      author='Ptdtan',
      author_email='ptdtan@gmail.com',
      license='MIT',
      packages=['zeus'],

      install_requires=['pandas','argparse','sklearn', 'scipy', 'statsmodels'],
      scripts=['src/zeus.py','src/logger.py'],
      zip_safe=False)
