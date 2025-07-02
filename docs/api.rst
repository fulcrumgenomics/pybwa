============
Introduction
============

Pybwa is a python module that makes it easy to align sequence
data.  It is a lightweight wrapper of `bwa`_.

This page provides a quick introduction in using `pybwa`_ followed by the
API.

Examples
========

Two alignment commands are supported: :code:`bwa aln` with the :class:`~pybwa.BwaAln` object and :code:`bwa mem` with :class:`~pybwa.BwaMem` object.
The constructor of both objects require either (1) a path prefix of the index (typically the FASTA), or (2) an already
created :class:`~pybwa.BwaIndex` object.

.. code-block:: python

   from pybwa import BwaAln
   aln = BwaAln(prefix="/path/to/genome.fasta")

or

.. code-block:: python

   from pybwa import BwaMem
   mem = BwaMem(prefix="/path/to/genome.fasta")

The :class:`~pybwa.BwaIndex` object is useful when re-using the same index, such that it only needs to be loaded into
memory once.  Both constructors for the :class:`~pybwa.BwaAln` and :class:`~pybwa.BwaMem` objects accept an index.

The :meth:`pybwa.BwaAln.align` method accepts a list of reads (as either strings or :class:`pysam.FastxRecord` s) to
align and return a *single* :class:`pysam.AlignedSegment` per input read:

.. code-block:: python

   for rec in aln.align(queries=["GATTACA"]):
      print(rec)

gives:

.. code-block:: console

   read.1	0	chr1	1	37	7M	*	0	0	GATTACA	*	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:7

The :meth:`pybwa.BwaMem.align` method accepts a list of reads (as either strings or :class:`pysam.FastxRecord` s) to
align and return a *list* of :class:`pysam.AlignedSegment` per input read:

.. code-block:: python

   for recs in mem.align(queries=["CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA"]):
      for rec in recs:
         print(rec)

gives:

.. code-block:: console

   chr1	0	chr1	1	60	35M	*	0	0	CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA	*	NM:i:0	MD:Z:35	AS:i:35	XS:i:0

The :meth:`pybwa.BwaAln.align` method accepts custom options provided as a :class:`~pybwa.BwaAlnOptions` object.
It is constructed directly with options set on the object:

.. code-block:: python

   opt = BwaAlnOptions()
   opt.max_mismatches = 5
   recs = aln.align(queries=["GATTACA"], opt=opt)


Similarly, the :meth:`pybwa.BwaMem.align` method accepts custom options provided as a :class:`~pybwa.BwaMemOptions` object.
It is constructed directly with options set on the object:

.. code-block:: python

   opt = BwaMemOptions()
   opt.min_seed_len = 32
   recs = aln.align(queries=["GATTACA"], opt=opt)

Note: the :meth:`~pybwa.BwaMemOptions.finalize` method will both apply the presets as specified by the
:meth:`~pybwa.BwaMemOptions.mode` option, as well as scale various other options (:code:`-TdBOELU`) based on the
:attr:`~pybwa.BwaMemOptions.match_score`.  The presets and scaling will only be applied to other options that have not
been modified from their defaults.  After calling the :meth:`~pybwa.BwaMemOptions.finalize` method, the options are
immutable, unless :code:`copy=True` is passed to :meth:`~pybwa.BwaMemOptions.finalize` method, in which case a copy
of the options are returned by the method.   Regardless, the :meth:`~pybwa.BwaMemOptions.finalize` method *does not*
need to be called before the :meth:`pybwa.BwaMem.align` is invoked, as the latter will do so (making a local copy).

API versus Command-line Differences
===================================

The reported alignments from `pybwa` may differ from those reported by the `bwa` command line.
This is true when `bwa` is run with a different number of threads (see :code:`bwa aln -t` and :code:`bwa mem -t`),
or when using a different number of chunks (see :code:`bwa mem -K`).

Finally, the following additions have been made to :code:`bwa aln/samse`:

#. The standard SAM tag :code:`HN` is added.  This is useful if we find too many hits
   (:attr:`~pybwa.BwaAlnOptions.max_hits`) and therefore no hits are reported in the :code:`XA` tag, we can still
   know how many were found.
#. The :py:attr:`~pybwa.BwaAlnOptions.with_md` option will add the standard SAM tag :code:`MD` to the :code:`XA` tag, 
   otherwise :code:`.` will be used.  This provides additional information on the quality of alternative alignments.

===
API
===

Bwa Aln
=======

.. autoclass:: pybwa.BwaAlnOptions
   :members:

.. autoclass:: pybwa.BwaAln
   :members:

.. autoclass:: pybwa.AuxHit
   :members:

.. autofunction:: pybwa.to_xa_hits

Bwa Mem
=======

.. autoclass:: pybwa.BwaMemOptions
   :members:

.. autoenum:: pybwa.BwaMemMode
   :members:

.. autoclass:: pybwa.BwaMem
   :members:

Bwa Index
=========

.. autoclass:: pybwa.BwaIndex
   :members:

.. autoenum:: pybwa.BwaIndexBuildMethod
   :members:
