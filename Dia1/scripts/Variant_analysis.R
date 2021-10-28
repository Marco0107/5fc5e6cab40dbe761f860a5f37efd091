library(data.table)
library(dplyr)
library(tibble)
setwd('~/ligandoma/LBB/')

# Importar as variantes preditas
vcf <- fread(cmd = "grep -v '##' results/amostra-lbb.HaplotypeCaller.filtered.vcf")
# Criar referencia genomica unica para comparar as variantes
vcf = vcf %>% mutate(COORD=paste(`#CHROM`, POS, sep='_'))

# Filtrar variantes que passaram no filtro do HaplotypeCaller
vcf_pass = dplyr::filter(vcf, FILTER=='PASS')
# Variantes filtradas / Total de variantes
print(paste0('Variantes filtradas: ', nrow(vcf_pass), '; ',
      paste0('Total de variantes: ', nrow(vcf)), '; ',
      paste0('(', round(nrow(vcf_pass) / nrow(vcf), digits = 4)*100, '%)') ))

# Importar as variantes referencia, padrao ouro
gabarito <- fread('data/pequeno-gabarito.vcf')
# Criar referencia genomica unica para comparar as variantes
gabarito = gabarito %>% mutate(COORD=paste(`#CHROM`, POS, sep='_'))

# Checar se as variantes do gabarito estao entre as variantes detectadas
all(gabarito$COORD %in% vcf_pass$COORD)
# Confirmar que o genotipo bate com o do padrao ouro
gt = vcf_pass[vcf_pass$COORD %in% gabarito$COORD,] %>%
  pull(`amostra-lbb`) %>% gsub('(\\d\\/\\d):\\S+', '\\1', .)
gt
gabarito$`AMOSTRA-LBB`

# Salvar o vcf filtrado:
write.table(vcf_pass, 'results/amostra-lbb.HaplotypeCaller.PASS.vcf', row.names = F, quote = F, sep = '\t')
