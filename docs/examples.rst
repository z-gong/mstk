Examples
========

Assign atom types using customized typing rule
----------------------------------------------

Atom type is the most fundamental element in a force field.
It determines which parameters should be used for describing the potential energy of each topological element (bond, angle, etc...) in a simulation system.
The assignment of atom types is usually the most cumbersome task in a simulation workflow.
For bio simulations, template matching is commonly used, because almost all bio-polymers are made of several kinds of repeating units.
However, the template matching won't work well for other molecules because there is no common repeating units.

`mstk` provide a typing engine :class:`~mstk.forcefield.typer.ZftTyper` for assigning atom types.
It is based on local chemical environment defined by SMARTS pattern and a hierarchical rule described in a type definition file.

In order to write the type definition file correctly, knowledge about
`SMARTS pattern <https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html>`_ is required.
The idea of the hierarchical rule is described in
`this article <https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.24244>`_.

Herein, an example of using :class:`~mstk.forcefield.typer.ZftTyper` to assign atom types is provided.
A type definition file should be fed to this typing engine.
The text below shows the type definition for linear alkane in `TEAM` force field.
Let's copy its content and save it to file `alkane.zft`.

.. literalinclude:: examples/alkane.zft

The script below will assign atom types for butane using this type definition.

.. literalinclude:: examples/type.py

A predefined typer :class:`~mstk.forcefield.typer.typer_primitive` is shipped with `mstk` for demonstration purpose.

Generate input files for simulation engines
-------------------------------------------

Note: this example requires that Packmol is installed.

In order to generate input files for simulation engines, a force filed file is required.
Currently, four formats are supported for storing force field:

1. `.zfp` - XML format used by `mstk`. Refer to :class:`~mstk.forcefield.Zfp` for details.
2. `.zff` - user-friendly text format used by `mstk`. Refer to :class:`~mstk.forcefield.Zff` for details.
3. `.ppf` - format used by DFF program. Refer to :class:`~mstk.forcefield.Ppf` for details.
4. `.ff` - format used by Agilio Padua's `fftool`. Refer to :class:`~mstk.forcefield.Padua` for details.

Herein, an example force field file for linear alkanes is given in `ZFF` format.
The parameters are taken from `TEAM_MGI` force field.
Let's copy its content and save it to file `alkane.zff`.

.. literalinclude:: examples/alkane.zff

The equivalent force field in `ZFP` format is

.. literalinclude:: examples/alkane.zfp

The workflow is as follows:

.. literalinclude:: examples/export.py

