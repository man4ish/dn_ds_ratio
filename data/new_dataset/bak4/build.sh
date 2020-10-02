cp sample.fa snpEff/data/kbase_v1/sequences.fa
cp sample_file.gff snpEff/data/kbase_v1/genes.gff
cd snpEff
java -jar snpEff.jar build -gff3 -v kbase_v1
java -Xmx8g -jar snpEff.jar kbase_v1 ../sample.vcf > ../sample.ann.vcf
cd ../../..
python parse_data.py
