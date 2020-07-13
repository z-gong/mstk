General applications
====================

Assign atom types using typing engine
-------------------------------------

`mstools` provide two typing engines for assigning atom types.
The first one :class:`~mstools.forcefield.typer.ZftTyper` is based on local chemical environment
described by SMARTS pattern and a hierarchical rule.
A type definition should be fed to this typing engine, either by file name or by string of its content.
In order to construct the type definition file correctly, knowledge about SMARTS is required.
https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html.
The idea of the hierarchical rule is described in http://dx.doi.org/10.1002/jcc.24244.
Here is an example of assigning OPLS atom types for butane and benzene.

.. literalinclude:: type.py

Generate input files for simulation engines
-------------------------------------------

Note: this example requires that Packmol is installed.

In order to generate input files for simulation engines, a force filed file is required.
This is an example force field file for linear alkanes in :class:`~mstools.forcefield.Zfp` format. The parameters are taken from `TEAM_MGI` force field.
Let's copy its content and save it to file `alkane.zfp`.

.. literalinclude:: alkane.zfp

The workflow is as follows:

.. literalinclude:: export.py
