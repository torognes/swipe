# bash completion for swipe        -*- shell-script -*-
#
# Provides tab-completion of swipe's options, of their values when they
# are names or small sets of numbers, and of file arguments.
#
# Installation (handled automatically by `make install`):
#   copy this file to one of the bash-completion lookup directories, e.g.
#     /usr/share/bash-completion/completions/swipe
#   or source it from your ~/.bashrc:
#     source /path/to/swipe.bash

_swipe()
{
    local cur prev
    COMPREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # Fallback for systems without the bash-completion package: a minimal
    # _filedir that completes files and directories.
    if ! declare -F _filedir >/dev/null 2>&1; then
        _filedir()
        {
            COMPREPLY=( $(compgen -f -- "$cur") )
        }
    fi

    # All recognised options, short and long.
    local all_opts="\
-a -b -c -C -d -D -e -E -F -G -h -H -i -I -k -K -m -M -N -o -p -q -Q -r \
-S -u -v -x -z \
--db --query --matrix --penalty --reward --gapopen --gapextend --strand \
--num_descriptions --num_alignments --min_score --max_score --evalue \
--minevalue --num_threads --outfmt --symtype --taxidlist --taxid \
--comp_based_stats --query_gencode --db_gencode --filter --subalignments \
--dump --out --dbsize --show_gis --show_taxid --help --version"

    case "$prev" in
        # Database: the base names of the BLAST databases (index files and
        # alias files, without their extension).
        -d|--db)
            local base
            COMPREPLY=( $(for base in $(compgen -f -- "$cur"); do
                              case "$base" in
                                  *.pin|*.nin|*.pal|*.nal) echo "${base%.*}" ;;
                              esac
                          done | sort -u) )
            COMPREPLY+=( $(compgen -d -- "$cur") )
            return 0
            ;;
        # Options whose argument is a file name -> complete with files.
        -i|--query|\
        -o|--out|\
        -x|--taxidlist|--taxid)
            _filedir
            return 0
            ;;
        # Score matrix: a built-in name, or a file name.
        -M|--matrix)
            COMPREPLY=( $(compgen -W "BLOSUM45 BLOSUM50 BLOSUM62 BLOSUM80 \
BLOSUM90 PAM30 PAM70 PAM250" -- "$cur") )
            local files=( $(compgen -f -- "$cur") )
            COMPREPLY+=( "${files[@]}" )
            return 0
            ;;
        -p|--symtype)
            COMPREPLY=( $(compgen -W "blastn blastp blastx tblastn tblastx \
sound" -- "$cur") )
            return 0
            ;;
        -S|--strand)
            COMPREPLY=( $(compgen -W "plus minus both" -- "$cur") )
            return 0
            ;;
        -m|--outfmt)
            COMPREPLY=( $(compgen -W "0 7 8 9 99" -- "$cur") )
            return 0
            ;;
        -N|--dump)
            COMPREPLY=( $(compgen -W "0 1 2" -- "$cur") )
            return 0
            ;;
        -Q|--query_gencode|\
        -D|--db_gencode)
            COMPREPLY=( $(compgen -W "1 2 3 4 5 6 9 10 11 12 13 14 15 16 \
21 22 23" -- "$cur") )
            return 0
            ;;
        -C|--comp_based_stats|\
        -F|--filter)
            COMPREPLY=( $(compgen -W "F" -- "$cur") )
            return 0
            ;;
        # Options whose argument is a number -> no value suggestion.
        -a|--num_threads|\
        -b|--num_alignments|\
        -c|--min_score|\
        -e|--evalue|\
        -E|--gapextend|\
        -G|--gapopen|\
        -k|--minevalue|\
        -K|--subalignments|\
        -q|--penalty|\
        -r|--reward|\
        -u|--max_score|\
        -v|--num_descriptions|\
        -z|--dbsize)
            return 0
            ;;
    esac

    # Complete option names (swipe has no positional argument).
    COMPREPLY=( $(compgen -W "$all_opts" -- "$cur") )
    return 0
}
complete -F _swipe swipe
