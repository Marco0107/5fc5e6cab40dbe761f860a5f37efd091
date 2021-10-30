**Importando vcf anotado com o Variant Effect Predictor (VEP). Esse
programa permite anotar frequências alélicas de variantes em diferentes
populações, anotar a consequência da variante para a proteína, dentre
outros.**

    vcf <- fread('../results/amostra-lbb.HaplotypeCaller.PASS.InBED.Otima.annotated.vcf', skip = '##')

### Obtenha a razão Ti/Tv (transitions e transversions) das variantes encontrada no cromossomo 22.

**Identificando e filtrando as variantes de base única (SNV)**

    snv <- paste0(vcf$REF, vcf$ALT)
    ind = sapply(snv, nchar) == 2 # indice para manter somente SNV (1+1 nt)
    snv = snv[ind]
    vcf_Ti_Tv = vcf[ind,] # para salvar e enviar o vcf
    write.table(vcf_Ti_Tv, 'amostra-lbb_Ti_Tv.vcf', quote = FALSE, row.names = FALSE, sep='\t')

**Somando todas as Transitions (Ti)**  
Transitions são mudanças: A &lt;-&gt; G ou C &lt;-&gt; T  
Vamos somar essas ocorrências

    Ti <- grep('AG|GA', snv) %>% length() + grep('CT|TC', snv) %>% length()

**Somando todas as Transversions (Tv)**  
Transversions são mudanças: A &lt;-&gt; C ou G &lt;-&gt; T ou G
&lt;-&gt; C ou A &lt;-&gt; T  
Soma-se todos os casos possíveis.

    Tv <- grep('AC|CA', snv) %>% length() + grep('GT|TG', snv) %>% length()
    Tv = Tv + grep('GC|CG', snv) %>% length() + grep('AT|TA', snv) %>% length()

**Razão Ti/Tv:**  
3.1666667

------------------------------------------------------------------------

### Quantas variantes são encontradas na região de 16000000 a 20000000?

Filtrando variantes que estejam contidas nessa região usando ocampo
‘POS’ do vcf.

    vcf_sub <- vcf %>% dplyr::filter(POS > 16000000 & POS < 20000000)

**Variantes na região:**  
59

------------------------------------------------------------------------

### Exiba o conteúdo da linha do VCF relativa a uma variante:

-   Não-sinônima..

As variantes não sinônimas são anotadas como ‘missense’ no arquivo vcf.
Por vezes elas não serão as únicas consequências para uma variante, pois
uma variante pode causar uma variante não sinônima em uma isoforma da
proteína e uma variante em região intronica se considerado outra
isoforma.

    vcf %>% dplyr::filter(grepl('missense', INFO)) %>% head(1)

    ##    #CHROM      POS ID REF ALT    QUAL FILTER
    ## 1:  chr22 16591593  .   A   G 4117.03   PASS
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             INFO
    ## 1: AC=2;AF=1.00;AN=2;DP=127;ExcessHet=3.0103;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=40.73;QD=32.42;SOR=0.916;CSQ=G|missense_variant|MODERATE|CCT8L2|ENSG00000198445|Transcript|ENST00000359963|protein_coding|1/1||ENST00000359963.4:c.958T>C|ENSP00000353048.3:p.Trp320Arg|1218|958|320|W/R|Tgg/Cgg|rs2236639||-1||SNV|HGNC|HGNC:15553|YES||P1|CCDS13738.1|ENSP00000353048|Q96SF2||UPI000006CF87||tolerated(1)|benign(0)|hmmpanther:PTHR11353&hmmpanther:PTHR11353:SF100&Pfam_domain:PF00118&Gene3D:3.50.7.10&Superfamily_domains:SSF52029|||0.8704|0.9501|0.8386|0.6875|0.9344|0.908|0.9421|0.9323|0.8835|0.9462|0.7833|0.9542|0.6521|0.8699|0.928|0.9082|0.9128|0.9542|gnomAD_ASJ||||20887961||||,G|downstream_gene_variant|MODIFIER|FABP5P11|ENSG00000240122|Transcript|ENST00000430910|processed_pseudogene||||||||||rs2236639|3496|-1||SNV|HGNC|HGNC:19328|YES||||||||||||||0.8704|0.9501|0.8386|0.6875|0.9344|0.908|0.9421|0.9323|0.8835|0.9462|0.7833|0.9542|0.6521|0.8699|0.928|0.9082|0.9128|0.9542|gnomAD_ASJ||||20887961||||
    ##            FORMAT                 amostra-lbb          COORD  DP FS    MQ    QD
    ## 1: GT:AD:DP:GQ:PL 1/1:0,127:127:99:4131,380,0 chr22_16591593 127  0 40.73 32.42

-   Variante no gnomAD v3.1.1 com MAF &lt; 0.01.

**Separar cada consequencia na variante por linha, uma variante pode
gerar multiplas anotações/consequencias na proteina**  
Isso é importante para capturar o campo correto das frequências
alélicas.

    vcf_separated <- vcf %>% 
      mutate(INFO_METRICS=gsub('(\\S+CSQ=[ATCG\\-]+)\\|\\S+','\\1', INFO),
             INFO=gsub('\\S+CSQ=([ATCG\\-]+\\|\\S+)','\\1', INFO))
    vcf_separated <- vcf_separated %>% separate_rows(INFO, sep = ',')

**Encontrando o campo da frequencia alelica (AF) no gnomAD**  
Aqui foi preciso identificar o campo que contém a informação de
interesse dentre os diferentes campos anotados pelo VEP.

    format = c('Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE')
    format = unlist(strsplit(format, split = '\\|'))
    ind=grep('gnomAD_AF$', format) # campo: 47

    vcf_separated = vcf_separated %>%
      mutate(gnomAD_AF=sapply(strsplit(INFO, split = '\\|'), '[', ind))

**Imprimindo um dos campos com MAF &lt; 0.01 no gnomAD**

    vcf_separated %>% 
      dplyr::filter(gnomAD_AF<0.01 & gnomAD_AF >=0) %>%
      dplyr::filter(grepl('missense', INFO)) %>% # opcional meu
      head(1) %>% as.data.frame()

    ##   #CHROM      POS ID REF ALT   QUAL FILTER
    ## 1  chr22 20426187  .   C   T 1749.6   PASS
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                INFO
    ## 1 T|missense_variant|MODERATE|SCARF2|ENSG00000244486|Transcript|ENST00000622235|protein_coding|11/11||ENST00000622235.4:c.1789G>A|ENSP00000477564.1:p.Glu597Lys|1860|1789|597|E/K|Gag/Aag|rs201929223||-1||SNV|HGNC|HGNC:19869||1|A2|CCDS46666.1|ENSP00000477564|Q96GP6||UPI000004D28D|1|tolerated(0.27)|benign(0.027)|mobidb-lite&hmmpanther:PTHR24043:SF5&hmmpanther:PTHR24043|||0.0034|0|0.0245|0|0|0|||0.002195|0|0.01182|0|0|0|0|0.002394|0|0.0245|AMR||||||||
    ##           FORMAT                  amostra-lbb          COORD  DP    FS MQ    QD
    ## 1 GT:AD:DP:GQ:PL 0/1:66,68:134:99:1757,0,1836 chr22_20426187 134 4.247 60 13.06
    ##                                                                                                                                                         INFO_METRICS
    ## 1 AC=1;AF=0.500;AN=2;BaseQRankSum=-2.643;DP=134;ExcessHet=3.0103;FS=4.247;MLEAC=1;MLEAF=0.500;MQ=60.00;MQRankSum=0.000;QD=13.06;ReadPosRankSum=0.548;SOR=0.420;CSQ=T
    ##   gnomAD_AF
    ## 1  0.002195
