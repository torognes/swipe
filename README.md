SWIPE
=====

Smith-Waterman database searches with inter-sequence SIMD parallelisation

SWIPE is a tool for performing rapid local alignment searches in amino
acid or nucleotide sequence databases. It is a highly optimized
implementation of the Smith-Waterman algorithm, using the SIMD parallel
computing instructions of common x86-64 CPUs (SSE2 is required, SSSE3
is faster): the optimal local alignment score is computed for every
database sequence, with no heuristic.

The method and a performance evaluation are described in detail in the
following publication:

Rognes T (2011)  
Faster Smith-Waterman database searches with inter-sequence SIMD
parallelisation.  
BMC Bioinformatics 12, 221.
doi:[10.1186/1471-2105-12-221](https://doi.org/10.1186/1471-2105-12-221)


## Quick start

SWIPE searches BLAST databases, prepared with `makeblastdb` (NCBI
BLAST+). Databases in format version 4 and version 5 (the default of
current `makeblastdb` versions, and the format of NCBI's pre-formatted
databases) are accepted:

```sh
makeblastdb -in proteins.fasta -dbtype prot -out proteins
swipe -d proteins -i queries.fasta -a 8 -e 0.001
```

The search type is chosen with `-p`: blastn (`-p 0`, nucleotide
query and database), blastp (`-p 1`, the default), blastx (`-p 2`,
translated query), tblastn (`-p 3`, translated database), tblastx
(`-p 4`, both translated). The output format is chosen with `-m`: plain
text (`-m 0`, the default, similar to BLAST's), XML (`-m 7`), tabular
(`-m 8`, the 12 columns of BLAST's `-m 8`, or `-m 9` with comment
lines).

Use `swipe -h` to get a short help, or see the [manual
page](https://github.com/torognes/swipe/blob/master/man/swipe.1) (`man
./man/swipe.1` before installation) for a complete description of the
options, the default gap penalties of each score matrix, and the
output formats. The file
[scoring.pdf](https://github.com/torognes/swipe/blob/master/scoring.pdf)
lists the scoring systems for which E-values can be computed.


## Install

Pre-compiled binaries for Linux (x86-64) are available on the
[releases page](https://github.com/torognes/swipe/releases) (with the
manual page and the shell completions from version 2.2.0 on). SWIPE can
also be installed with
Homebrew:

```sh
brew install brewsci/bio/swipe
```

To compile SWIPE from source, you need a C++11 compiler (GCC 4.8.5 or
later, or clang) and GNU make:

```sh
git clone https://github.com/torognes/swipe.git
cd swipe/
make
# or, with clang
make CXX=clang++
```

`make install` copies the binary, the manual page and the bash and
zsh completions under `/usr/local` (`make install PREFIX=$HOME/.local`
to install elsewhere; `DESTDIR` is supported), and `make uninstall`
removes them. The usual variables `CXX`, `CXXFLAGS`, `CPPFLAGS` and
`LDFLAGS` are honoured (e.g. `make CXXFLAGS=-O2`).


## Shell auto-completion

`make install` installs the completions for bash
(`share/bash-completion/completions/swipe`) and zsh
(`share/zsh/site-functions/_swipe`). Without installation, they can be
loaded from the source tree:

```sh
# bash (add to ~/.bashrc)
source /path/to/swipe/completion/swipe.bash

# zsh (add the directory to your $fpath, before `compinit`, in ~/.zshrc)
fpath=(/path/to/swipe/completion $fpath)
```


## Versions

The source code of every version from 2.0.5 onward is available on
[GitHub](https://github.com/torognes/swipe/tags), and pre-compiled
binaries of versions 2.1.0 and later on the
[releases page](https://github.com/torognes/swipe/releases). Versions
1.0 to 2.0.4 are no longer distributed. The changes of each version are
listed in the
[CHANGES](https://github.com/torognes/swipe/blob/master/CHANGES) file,
and the [README](https://github.com/torognes/swipe/blob/master/README)
gives more details on the options and examples.

The MPI version of SWIPE (mpiswipe) was removed in version 2.2.0
(version 2.1.2 is the last one that includes it).


## License and contact

SWIPE is distributed under the GNU Affero General Public License,
version 3. Please report bugs and suggestions on the [issue
tracker](https://github.com/torognes/swipe/issues).
