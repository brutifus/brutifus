.. _faq:

FAQ 
===

**Q**: What's with `brutus <https://github.com/fpavogt/brutus>`_? Is it the same as brutifus?

   **A**: the brutus code is the ancestor to brutifus, and is now deprecated. It may be that 
   not all the functionalities of brutus will be implemented inside brutifus: it'll all 
   depend on what I use brutifus for.

**Q**: I'm getting some weird LaTeX error massage ... what is going ?

   **A**: I'm assuming brutifus is run on a machine with a proper 
   [`citation needed <https://en.wikipedia.org/wiki/Wikipedia:Citation_needed>`_] 
   system-wide LaTeX installation, to get prettier plots. Although I very much *do not 
   recommend it*, you can by-pass this problem by running brutifus with the ``--no-systemtex``
   flag::
   
       >>> run brutus_execute --params pickle_fn --no-systemtex