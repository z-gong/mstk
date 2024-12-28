import os
import setuptools

with open('README.md') as f:
    long_description = f.read()

dir_scripts = 'scripts'
scripts = [os.path.join(dir_scripts, f) for f in os.listdir(dir_scripts)]

setuptools.setup(
    name='mstk',
    version='0.3.22',
    author='Zheng Gong',
    author_email='z.gong@outlook.com',
    description='Molecular simulation toolkit',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/z-gong/mstk',
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={
        'mstk': ['data/forcefield/*']
    },
    scripts=scripts,
    python_requires='>=3.6',
    # dependency can be a mess if conda and pip are mixed
    # better let user install requirements by themselves
    install_requires=[],
    classifiers=[
        'License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)',
        'Operating System :: OS Independent'
    ]
)
