============
Introduction
============

Pybwa is a python module that makes it easy to align sequence
data.  It is a lightweight wrapper of `bwa`_.

This page provides task-oriented examples followed by the full API reference.

========
Examples
========

Align short reads (bwa aln)
===========================

Use :class:`~pybwa.BwaAln` for aligning short reads (<100bp, e.g. Illumina).  The
:meth:`~pybwa.BwaAln.align` method returns one :class:`~pysam.AlignedSegment` per input read.

.. code-block:: python

   from pybwa import BwaAln

   aln = BwaAln(prefix="/path/to/genome.fasta")
   for rec in aln.align(queries=["GATTACA"]):
       print(rec)

gives:

.. code-block:: console

   read.1	0	chr1	1	37	7M	*	0	0	GATTACA	*	XT:A:U	NM:i:0	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:7


Align reads (bwa mem)
=====================

Use :class:`~pybwa.BwaMem` for aligning longer reads (>70bp, Illumina 100bp+, PacBio, ONT).  The
:meth:`~pybwa.BwaMem.align` method returns a *list* of :class:`~pysam.AlignedSegment` per input
read, since a single query may produce primary, supplementary, and secondary alignments.

.. code-block:: python

   from pybwa import BwaMem

   mem = BwaMem(prefix="/path/to/genome.fasta")
   for recs in mem.align(queries=["CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA"]):
       for rec in recs:
           print(rec)

gives:

.. code-block:: console

   chr1	0	chr1	1	60	35M	*	0	0	CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA	*	NM:i:0	MD:Z:35	AS:i:35	XS:i:0


Customize alignment options
===========================

Both aligners accept an options object to control alignment behavior.

For :class:`~pybwa.BwaAln`, construct a :class:`~pybwa.BwaAlnOptions` and pass it to
:meth:`~pybwa.BwaAln.align`:

.. code-block:: python

   from pybwa import BwaAln, BwaAlnOptions

   aln = BwaAln(prefix="/path/to/genome.fasta")
   opt = BwaAlnOptions()
   opt.max_mismatches = 5
   recs = aln.align(queries=["GATTACA"], opt=opt)

For :class:`~pybwa.BwaMem`, construct a :class:`~pybwa.BwaMemOptions` and pass it to
:meth:`~pybwa.BwaMem.align`:

.. code-block:: python

   from pybwa import BwaMem, BwaMemOptions

   mem = BwaMem(prefix="/path/to/genome.fasta")
   opt = BwaMemOptions()
   opt.min_seed_len = 32
   recs = mem.align(queries=["CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA"], opt=opt)

Note: the :meth:`~pybwa.BwaMemOptions.finalize` method will both apply the presets as specified by the
:attr:`~pybwa.BwaMemOptions.mode` option, as well as scale various other options (:code:`-TdBOELU`) based on the
:attr:`~pybwa.BwaMemOptions.match_score`.  The presets and scaling will only be applied to other options that have not
been modified from their defaults.  After calling the :meth:`~pybwa.BwaMemOptions.finalize` method, the options are
immutable, unless :code:`copy=True` is passed to :meth:`~pybwa.BwaMemOptions.finalize` method, in which case a copy
of the options are returned by the method.  Regardless, the :meth:`~pybwa.BwaMemOptions.finalize` method *does not*
need to be called before the :meth:`pybwa.BwaMem.align` is invoked, as the latter will do so (making a local copy).


Reuse an index across multiple alignments
==========================================

Load a :class:`~pybwa.BwaIndex` once and pass it to multiple aligners to avoid re-loading
the index into memory for each aligner:

.. code-block:: python

   from pybwa import BwaIndex, BwaAln, BwaMem

   index = BwaIndex(prefix="/path/to/genome.fasta")
   aln = BwaAln(index=index)
   mem = BwaMem(index=index)

This is also useful when aligning multiple batches of reads; construct the aligner once and call
:meth:`~pybwa.BwaAln.align` or :meth:`~pybwa.BwaMem.align` repeatedly.


Parse alternative hits from the XA tag
=======================================

BWA reports alternative alignment locations in the :code:`XA` SAM tag.  Use
:func:`~pybwa.to_xa_hits` to parse these into :class:`~pybwa.AuxHit` objects:

.. code-block:: python

   from pybwa import BwaAln, to_xa_hits

   aln = BwaAln(prefix="/path/to/genome.fasta")
   rec = aln.align(queries=["GATTACA"])[0]
   hits = to_xa_hits(rec)
   for hit in hits:
       print(f"{hit.refname}:{hit.start}-{hit.end} edits={hit.edits}")

You can also parse an XA tag value directly (without the leading :code:`XA:Z:` prefix):

.. code-block:: python

   hits = to_xa_hits("chr4,-97592047,24M,3;chr8,-32368749,24M,3;")


Control BWA verbosity
=====================

By default BWA outputs informational and warning messages to stderr.  Use
:func:`~pybwa.set_bwa_verbosity` to suppress or increase output:

.. code-block:: python

   from pybwa import BwaVerbosity, set_bwa_verbosity

   set_bwa_verbosity(BwaVerbosity.QUIET)  # suppress all output


Use long-read presets (PacBio/ONT)
==================================

Use :class:`~pybwa.BwaMemMode` to apply presets for PacBio or Oxford Nanopore reads.  This adjusts
multiple options (penalties, seed length, etc.) to be appropriate for the given read type:

.. code-block:: python

   from pybwa import BwaMem, BwaMemOptions, BwaMemMode

   mem = BwaMem(prefix="/path/to/genome.fasta")
   opt = BwaMemOptions(mode=BwaMemMode.PACBIO)
   recs = mem.align(queries=["ACGT" * 100], opt=opt)

Available presets:

- :attr:`~pybwa.BwaMemMode.PACBIO` -- PacBio CLR reads
- :attr:`~pybwa.BwaMemMode.ONT2D` -- Oxford Nanopore 2D reads
- :attr:`~pybwa.BwaMemMode.INTRACTG` -- intra-species contig alignment


API versus Command-line Differences
====================================

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

Verbosity
=========

.. autoenum:: pybwa.BwaVerbosity
   :members:

.. autofunction:: pybwa.set_bwa_verbosity
