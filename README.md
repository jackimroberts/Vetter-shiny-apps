In this directory are R shiny apps to explore published data from the Vetter Lab.
Note that this does not include the data files but running the app will pull them from GEO.

h4(strong("Microglia during retinal development"))),
    strong("Whole retina timecourse, bulk RNA-seq:"),("Triplicate central retina samples at e16.5, and single whole retina samples at P7 and P60. 
                                                        Published in Anderson et al. Cell Reports, 2019. GSE123757"),
      tags$br(), #linebreak
    strong("Retinal microglia timecourse, bulk RNA-seq:"),("Three biological replicates of retinal microglia (GFP+, CD45+, CCR2âˆ’) at e16.5, P7 and P60. 
                                                             Published in Anderson et al. Cell Reports, 2019. GSE123757"),
      tags$br(),
    strong("Brain microglia timecoure, bulk RNA-seq:"),("Duplicate e16, P7 brain (CD45+CD11B+CX3CR1-GFP+) and P60 cortex microglia (CD45intCD11BintCX3CR1-GFP+) 
                                                          Fastq files from Matcovich Natan et al., Science, 2016 GSE79812. Also analyzed in Anderson et al. 
                                                          Cell Reports, 2019."),
      tags$br(),
    strong("Itgax+ microglia, bulk RNA-seq:"),("P5 CD11c-Cre-gfp;Rosa-TDTomato sorted for active Cd11c expression (GFP+TdTom+) or Cd11c lineage (GFP-TdTom+). This is unplished and won't show up"),
      tags$br(),
    strong("Axl MerTK KO, bulk RNA-seq:"),("Triplicate P4 retinal microglia from Mertk KO, Axl KO and P7 Mertk/Axl dKO.
                                        Published in Anderson et al. eLife, 2022. GSE192601"),
      tags$br(),
    strong("Retinal micorglia PLX and Bax cKO at P7, scRNA-seq:"),("FACS sorted from Bax WT and littermate KO (CD45+ CD11b+/GFP+ CCR2-), 
    PLX3397-treated CX3CR1-GFP/+ and Vehicle CX3CR1-GFP/+ (CD45+GFP+ Ly6C-). Published in Anderson et al. eLife, 2022. GSE192600")
tags$a(href="https://vetterdata.chpc.utah.edu/Zhang", h4(strong("Epigenetic regulation of retinal development"))),
    strong("Ezh2 cKO retina at E16.5, bulk RNA-seq:"), #bold
      ("Peripheral one-third of E16.5 retinas was dissected from four biological replicates of Ezh2"),tags$sup("fl/−"), #superscript
      ("CKO (Pax6-αCre) animals and Ezh2"),tags$sup("fl/+"),("controls from the same litter. Published in Zhang et al. Developmental Biology, 2015. GSE65082"),
      tags$br(),
    strong("Jarid2 cKO retina at E16.5, bulk RNA-seq:"),em("add details"),
      tags$br(),
    strong("Jarid2 cKO retina at P0, bulk RNA-seq:"),em("add details"),
      tags$br(),
    strong("Jarid2 cKO retina at e18.5, scRNA-seq:"),em("add details"),
      tags$br(),
    strong("Jarid2 cKO cortex at E15.5, bulk RNA-seq:"),("Singlet cortical hemispheres from Jarid2 fl/fl, Emx1"), tags$sup("cre/+"),
            (" (cKO) and Jarid2+/+, Emx1"),tags$sup("cre/+"),(" (control)"),
