# this file provides bash completors for some of binaries shipped with
# eWoms. This allows to implement tab completion of parameters for
# these binaries without having to replace the default completor

# find the directory of the ewoms_product_completors.sh file
FOO_DIR="."
if test -n "$BASH_SOURCE"; then
    FOO_DIR=$(dirname "$BASH_SOURCE") 
fi

source "$FOO_DIR/ewoms_libbash_completors.sh"
complete -o nospace -F _ewoms_parameter_completor ebos

unset FOO_DIR
