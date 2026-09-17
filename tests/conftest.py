"""Every test drives a compiled Stan program; the first use of each compiles
it (about ten seconds) and the module cache in ``dgsn.stan.run`` keeps it for
the session.  Run under the ``dgsn`` environment, where CmdStan is found."""
