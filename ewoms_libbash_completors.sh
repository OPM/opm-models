# this snippet enables flow-parameter completion via the tabulator key
# for bash. it can be used by adding the following
#
#   source $PATH_TO_EWOMS/ewoms_libbash_completors.sh
#   complete -D -o nospace -F _ewoms_generic_parameter_completor
#
# to the ~/.bashrc file. (Actually this works for any eWoms based
# simulator, not just for flow!) Note that, depending on how your
# operating system implements bash parameter completion, you might
# need to set the DEFAULT_COMPLETION_LOADER environment variable. The
# above has been verfied to work with Debian based linux distributions
# (i.e., Ubuntu >= 16.04) and openSUSE.


# this is a bash readline completer for the case where a binary is
# already known to be an eWoms simulator.
_ewoms_parameter_completor()
{
    if test "$COMP_WORDS" == ""; then
        return 0
    fi

    cmd="${COMP_WORDS[0]}"
    cur="${COMP_WORDS[COMP_CWORD]}"
    fullcmd="$(which "$cmd")"
    ALL_OPTS=$("$fullcmd" --help 2> /dev/null | grep '^ *--' | sed 's/ *\(--[a-zA-Z0-9\-]*\)=.*/\1=/')
    ALL_OPTS=$(echo "$ALL_OPTS" | sed 's/^ *--help.*/--help/')
    COMPREPLY=( $(compgen -A file -W "$ALL_OPTS" -- "${cur}") )
}

# this is a bash readline default completer which attempts to find out
# if a given binary is an eWoms simulation. this needs to be set as a
# default completer because the name of eWoms binaries cannot be known
# a-priori.
_ewoms_generic_parameter_completor()
{
    if test "$COMP_WORDS" == ""; then
        return 0
    fi

    COMPREPLY=()
    local cmd cur ALL_OPTS
    cmd="${COMP_WORDS[0]}"
    cur="${COMP_WORDS[COMP_CWORD]}"

    fullcmd="$(which "$cmd")"
    if test -z "$fullcmd" || \
       ! test -x "$fullcmd" || \
       (! test -f "$fullcmd" && ! test -h "$fullcmd" ) || \
       ! test -r "$fullcmd" || \
       ! grep -q "Ewoms[a-zA-Z0-9]*Simulator[a-zA-Z0-0]" "$fullcmd"
    then
        if test -n "$DEFAULT_COMPLETION_LOADER"; then
            "$DEFAULT_COMPLETION_LOADER" $@
        elif type -t _completion_loader 2>&1 > /dev/null; then
            # the default DEFAULT_COMPLETION_LOADER variable has not
            # been set and the _completion_loader function exists, so
            # we use _completion_loader as the default completer.
            _completion_loader $@
        else
            return 1
        fi

        return $?
    fi

    _ewoms_parameter_completor $@
    return 0
}
