#Fasta sequence to uppercase.
for f in *.fa;do
awk '/^>/ {print($0)}; /^[^>]/ {print(toupper($0))}' ${f} >tmp
mv tmp ${f}
done
