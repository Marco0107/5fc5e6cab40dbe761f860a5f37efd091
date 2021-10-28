Ordem de execução dos scripts:
1) quality.sh
2) gatk4.slurm
3) HaplotypeCaller.slurm
4) Variant_analysis.R


Observações e conclusões após cada etapa:
1) A qualidade das leituras é satisfatória e não for necessário nenhum processo de limpeza/aparamento.
2) Alinhamento ao genoma de referência contendo o chr22 conforme as boas práticas do GATK.
3) Identificação de variantes germline usando a ferramenta do GATK e posterior filtragem por parâmetros de qualidade e cobertura conforme: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering
4) Comparação com variantes fornecidas e filtragem das variantes que não passaram os filtros mínimos do GATK.
