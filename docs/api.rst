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
It is constructed directly and options set on the object:

.. code-block:: python

   opt = BwaAlnOptions()
   opt.max_mismatches = 5
   recs = aln.align(queries=["GATTACA"], opt=opt)


The :meth:`pybwa.BwaMem.align` method accepts custom options provided as a :class:`~pybwa.BwaMemOptions` object.
It is constructed via the :class:`~pybwa.BwaMemOptionsBuilder` class, to support scaling gap open and extend penalties
when a using custom match score, or the specification of presets (via `mode`).

.. code-block:: python

   builder = BwaMemOptionsBuilder()
   builder.min_seed_len = 32
   opt: BwaMemOptions = builder.build()
   recs = aln.align(queries=["GATTACA"], opt=opt)

The :class:`~pybwa.BwaIndex` object is useful when re-using the same index, such that it only needs to be loaded into memory
once.
Both constructors for the :class:`~pybwa.BwaAln` and :class:`~pybwa.BwaMem` objects accept an index.

API versus Command-line Differences
===================================

The reported alignments from `pybwa` may differ from those reported by the `bwa` command line.
In particular when the latter is run with multiple threads (see :code:`bwa aln -t` and :code:`bwa mem -t`),
and when the latter processes multiple chunks (see :code:`bwa aln -m` and :code:`bwa mem -K`).

Furthermore, when multiple alignments with the same *alignment score* exist for one read, both the chosen primary
alignment *AND* number of non-primary (i.e. secondary) alignments may differ.
As a result, the mapping quality may also differ slightly.
This is due an implementation detail in which the order of random numbers are used.

===
API
===

Bwa Aln
=======

.. autoclass:: pybwa.libbwaaln.BwaAlnOptions
   :members:

.. autoclass:: pybwa.libbwaaln.BwaAln
   :members:

Bwa Mem
=======

.. autoclass:: pybwa.BwaMemOptionsBuilder
   :members:

.. autoclass:: pybwa.BwaMemOptions
   :members:

.. autoclass:: pybwa.BwaMemMode
   :members:

.. autoclass:: pybwa.BwaMem
   :members:
