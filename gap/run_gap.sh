#!/bin/zsh

# -- create output directory ------------------------------------------

DIR="./gapout_triangle"
if ! [ -d "$DIR" ]; then
  mkdir $DIR 
fi

# -- store local path -------------------------------------------------

LOCAL_DIR=$(pwd)

# -- run gap ----------------------------------------------------------

alias gap=/home/fox/Git/gap-4.12.2/bin/gap.sh


gap -b -q << EOI

ChangeDirectoryCurrent("$LOCAL_DIR");
Read("./lins_triangle_group.g");

quit;

EOI
