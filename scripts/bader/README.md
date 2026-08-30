# Bundled Bader charge-analysis program

`bader.x` is the **Grid-Based Bader Analysis** program by the Henkelman group
(https://theory.cm.utexas.edu/henkelman/code/bader/), Version 1.05 (2023-08-19).

It is bundled with AbacusCopilot and invoked **transparently** by task **1005
(Bader Charge)** to partition the total charge density into atomic basins and
integrate per-atom Bader charges.  Users do not need to install or invoke
`bader.x` themselves.

## License

The Bader program is distributed under the GNU General Public License, which is
compatible with AbacusCopilot's GPL-3.0.  See the Henkelman group page for the
canonical terms.  Source is in `bader/`.

## Rebuilding (only needed if a platform's prebuilt binary is missing)

Requires `gfortran`.  From the `abacuscopilot/scripts/bader` directory:

```bash
# Linux: static link works
( cd bader && make -f makefile.lnx_gfortran ) && cp bader/bader bader.x

# macOS: static link is unavailable, drop it
( cd bader && make -f makefile.osx_gfortran LINK="" ) && cp bader/bader bader.x
```

`setup.sh` performs this automatically at install time (tries static first,
falls back to non-static).
