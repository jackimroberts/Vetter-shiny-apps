# Vetter Lab data

In this directory are R shiny apps to explore published data from the Vetter Lab.

All data (including unpublished) is on the secure shiny server, but I thought I'd offer these examples for people to play with or add new plots to.

## Microglia during retinal development
``` r
runApp('ShinyApp/Anderson')
```
**Whole retina timecourse, bulk RNA-seq**  
Triplicate central retina samples at e16.5, and single whole retina samples at P7 and P60.  
Published in [Anderson et al. Cell Reports, 2019](%2210.1016/j.celrep.2019.04.062). In GEO repository: [GSE123757](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123757)

**Retinal microglia timecourse, bulk RNA-seq**  
Three biological replicates of retinal microglia (GFP+, CD45+, CCR2-) at e16.5, P7 and P60.  
Published in Anderson et al. Cell Reports, 2019. In GEO repository: [GSE123757](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123757)

**Brain microglia timecourse, bulk RNA-seq**  
Duplicate e16, P7 brain (CD45+CD11B+CX3CR1-GFP+) and P60 cortex microglia (CD45-int CD11B-int CX3CR1-GFP+)
Fastq files from Matcovich Natan et al., Science, 2016 GSE79812. Also analyzed in Anderson et al. Cell Reports, 2019.

**Axl MerTK KO, bulk RNA-seq**  
Triplicate P4 retinal microglia from Mertk KO, Axl KO and P7 Mertk/Axl dKO.  
Published in [Anderson et al. eLife, 2022](https://doi.org/10.7554/elife.76564). In GEO repository: [GSE192601](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192601)

**Retinal micorglia PLX and Bax cKO at P7, scRNA-seq: FACS sorted from Bax**  
WT and littermate KO (CD45+ CD11b+/GFP+ CCR2-), PLX3397-treated CX3CR1-GFP/+ and Vehicle CX3CR1-GFP/+ (CD45+GFP+ Ly6C-). 
Published in Anderson et al. eLife, 2022. [GSE192600](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192600)

## Epigenetic regulation of retinal development
```r
runApp('ShinyApp/Zhang')
```
**Ezh2 cKO retina at E16.5, bulk RNA-seq**  
Peripheral one-third of E16.5 retinas was dissected from four biological replicates of Ezh2 fl/−, CKO (Pax6-αCre) animals and Ezh2 fl/+ controls from the same litter. 
Published in [Zhang et al. Developmental Biology, 2015](https://doi.org/10.1016/j.ydbio.2015.05.010). [GSE65082](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65082)

**Jarid2 cKO retina at E16.5, bulk RNA-seq.** Unpublished

**Jarid2 cKO retina at P0, bulk RNA-seq.** Unpublished

**Jarid2 cKO retina at e18.5, scRNA-seq.**  
Published in [Zhang et al. Cell Reports, 2023](https://doi.org/10.1016/j.celrep.2023.112237). [GSE202734](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202734)

**Jarid2 cKO cortex at E15.5, bulk RNA-seq**  
Single cortical hemispheres from Jarid2 fl/fl, Emx1 cre/+ (cKO) and Jarid2+/+, Emx1 cre/+ (control).
Unpublished
