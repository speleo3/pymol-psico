'''
psico installation script
'''

from distutils.core import setup

try:
	from psico import __version__
except:
	print('Warning: could not import version')
	__version__ = '0.0'

distribution = setup(
	name='psico',
	version=__version__,
	packages=[
		'psico',
		],
	author='Thomas Holder and Steffen Schmidt',
	url='https://github.com/speleo3/pymol-psico',
	description='Pymol ScrIpts COllection',
	license='BSD-2-Clause',
	requires=['pymol'],
	)
