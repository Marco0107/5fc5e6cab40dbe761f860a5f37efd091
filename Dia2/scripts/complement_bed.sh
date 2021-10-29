#!/usr/bin/bash

# Contar bases no chr do genoma de referencia
bases=$(grep -v '>' grch38.chr22.fasta | wc -m)

# Gerar um arquivo com o numero de bases do genoma para gerar o .bed complementar
echo -e 'chr22\t'$bases > my_genome.bed

# Obter as regioes nao cobertas usando o bedtools
bedtools complement -i coverage.bed -g my_genome.bed > coverage_complement.bed
