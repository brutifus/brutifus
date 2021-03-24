.. include:: ./substitutions.rst
.. _faq:

FAQ
===

**Q**: What's with `brutus <https://github.com/fpavogt/brutus>`_? Is it the same as brutifus?

   **A**: the brutus code is the ancestor to brutifus, and is now deprecated. One should note that
   not all the brutus capabilities were ported to brutifus.

**Q**: I'm getting some weird LaTeX error message ... what is going on ?

   **A**: If you set ``systemtex: True`` in ``params_brutifus.yaml``, the code assumes that it is being
   run on a machine with a proper [`citation needed <https://en.wikipedia.org/wiki/Wikipedia:Citation_needed>`_]
   system-wide LaTeX installation, which should lead to prettier plots. If this fails for you, you
   can a) fix your LaTeX installation, or b) disable its use in favor of matplotlib's native one by
   setting ``systemtex: False``.
