General applications
====================

Assign atom types using typing engine
-------------------------------------

`mstools` provide two typing engines for assigning atom types.

1. :class:`~mstools.forcefield.typer.ZftTyper` is based on local chemical environment
described by SMARTS pattern and a hierarchical rule described in a type definition file.

2. :class:`~mstools.forcefield.typer.DffTyper` calls `DFF` program to assign atom types.
(The dependence of `DFF` will be removed in the future)

Herein, an example of using :class:`~mstools.forcefield.typer.ZftTyper` to assign atom types is provided.
A type definition should be fed to this typing engine, either by file name or by string of its content.
In order to construct the type definition file correctly, knowledge about
`SMARTS pattern <https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`_ is required.
The idea of the hierarchical rule is described in
`this article <https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.24244>`_.
The script below will assign atom types for butane using the convention of `TEAM` force field.

.. literalinclude:: type.py


Generate input files for simulation engines
-------------------------------------------

Note: this example requires that Packmol is installed.

In order to generate input files for simulation engines, a force filed file is required.
Currently, three formats are supported for storing force field:

1. `.zfp` - default format used by `mstools`. Refer to :class:`~mstools.forcefield.Zfp` for details.
2. `.ppf` - format used by DFF program. Refer to :class:`~mstools.forcefield.Ppf` for details.
3. `.ff` - format used by Agilio Padua's `fftool`. Refer to :class:`~mstools.forcefield.Padua` for details.

Herein, an example force field file for linear alkanes is given in `ZFP` format.
The parameters are taken from `TEAM_MGI` force field.
Let's copy its content and save it to file `alkane.zfp`.

.. literalinclude:: alkane.zfp

The workflow is as follows:

.. literalinclude:: export.py
