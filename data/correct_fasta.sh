zcat Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz | sed '13242549,16089942!d' > chr6.txt
sed '901328,$!d' chr6.txt | sed 's/AGTCCAGGTCCCGGAAAGGAGGCAGGATGCTTGGGACAGGAATCGGAGGCGTTGGGGTGT/AGTCCAGGTCCCGGAAAGGGAGGCAGGATGCTTGGGACAGGAATCGGAGGCGTTGGGGTGT/' > chr6_after.txt
echo "" > chr6_after_first.txt
grep .$ -o chr6_after.txt >> chr6_after_first.txt
sed -i '$d' chr6_after_first.txt
paste -d "" chr6_after_first.txt chr6_after.txt > temp.after.txt
cat temp.after.txt | head -n -1 | sed 's/.$//' > old.after.txt
tail -n1 temp.after.txt >> old.after.txt

head -n901327 chr6.txt > new.chr6.txt
cat old.after.txt >> new.chr6.txt

zcat Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz | head -n13242548 > corrected_Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
cat new.chr6.txt >> corrected_Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
zcat Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz | sed '16089943,$!d' >> corrected_Sus_scrofa.Sscrofa11.1.dna.toplevel.fa
