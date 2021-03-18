# to remove mutations with chromosome in the format chr*_*
# --> remove all mutations with chromosome contains "_"

file=$1
ofile=$2

cat "$file" | grep "^#" >> "$ofile"  && cat "$file" | grep -v "^#" | awk '!($1 ~ /[_]/)' >> "$ofile"

