#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# installing bioconductor
library(BiocManager)
options(repos = BiocManager::repositories())

# loading and installing all the required libs
library(shiny)
library(shinydashboard)

# loading bioconductor packages
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)

# igvShiny
# library(igvShiny)

# loading heatmap
load(file = "data/ht_with_names_for_shiny.RData")
htComplete = draw(htComplete)

# creating header
header = dashboardHeader(title = "Halo Atlas")

# creating sidebar
sidebar = dashboardSidebar(
    sidebarMenu(
        menuItem("Interactive heatmap", tabName = "interactiveHeatMap"),
        menuItem("Annotation legend", tabName = "heatmapLegend"),
        menuItem("Genome browser", tabName = "igv"),
        menuItem("About & Help", tabName = "help")
    )
)

# creating body
body = dashboardBody(
    tabItems(
         tabItem(
            tabName = "interactiveHeatMap",
            fluidRow(
                box(
                    title = "Interactive heatmap",
                    width = 12, solidHeader = TRUE, status = "primary",
                    InteractiveComplexHeatmapOutput("ht_atlas",
                                                    layout = "1|2|3",
                                                    output_ui = NULL,
                                                    title3 = NULL,
                                                    width1 = 1000,
                                                    width2 = 1000,
                                                    height1 = 650,
                                                    height2 = 650)
                )
            )
        ),
        tabItem(
            tabName = "heatmapLegend",
            fluidRow(
                box(
                    title = "Annotation legend",
                    width = 12, height = 650, solidHeader = TRUE, status = "primary",
                    imageOutput("heatmapLegend")
                )
            )
        ),
        # tabItem(
        #     tabName = "igv",
        #     fluidRow(
        #         box(
        #             title = "Genome browser",
        #             width = 12, solidHeader = TRUE, status = "primary",
        #             igvShinyOutput("igvShiny_0")
        #         )
        #     )
        # ),
        
        # taken and modified from one of Paul Shannon's example of igv.js
        # https://github.com/alanlorenzetti/igvShiny/blob/master/igv.js.examples/arcs-bedpe.html
        # 
        tabItem(
            tabName = "igv",
            fluidRow(
                box(
                    title = "Genome browser (loading could take a few seconds)",
                    width = 12, solidHeader = TRUE, status = "primary",
                    HTML('
                    <div id="igvDiv" style="height: auto"></div>
                    
                    <script type="module">
                    
                    import igv from "https://cdn.jsdelivr.net/npm/igv@2.6.8/dist/igv.esm.js"

                    var options =
                    {
                    
                        locus: "NC_001869.1:17,638-27,205",
                    
                        reference: {
                            id: "Halobacterium salinarum NRC-1",
                            fastaURL: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/Hsalinarum.fa",
                            indexURL: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/Hsalinarum.fa.fai",
                            wholeGenomeView: false
                        },
                        
                        tracks: [
                            {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/ribosomal_RNA_TP2-fwd.bedgraph.gz",
                                type: "wig",
                                format: "bedgraph",
                                autoscaleGroup: "ribo",
                                name: "Ribo-Seq TP2 (+)",
                                color: "darkgrey"
                            },
                            {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/total-RNA-TP2-rev.bedgraph.gz",
                                type: "wig",
                                format: "bedgraph",
                                autoscaleGroup: "tot",
                                name: "RNA-Seq TP2 (+)",
                                color: "darkgrey"
                            },
                            
                            {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/TPS_gff_fwd.gff3",
                                type: "annotation",
                                format: "gff3",
                                name: "TPS (+)",
                                height: 30,
                                color: "#59A14F",
                                displayMode: "SQUISHED"
                            },
                            
                            {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/mergedBRs_interaction-regions-entire-genome-fwd.gff3",
                                type: "annotation",
                                format: "gff3",
                                name: "SmAP1 interaction (+)",
                                height: 30,
                                color: "#E15759",
                                displayMode: "SQUISHED"
                            },
                            
                         {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/Hsalinarum-846asRNAs-deAlmeida2019.gff3",
                                type: "annotation",
                                format: "gff3",
                                name: "Antisense RNAs",
                                color: "#B07AA1",
                                displayMode: "EXPANDED"
                            },
                        
                            {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/Hsalinarum-gene-annotation-pfeiffer2019-adjusted-names.gff3",
                                type: "annotation",
                                searchable: true,
                                format: "gff3",
                                name: "Genes",
                                color: "#4E79A7",
                                displayMode: "COLLAPSED"
                            },
                            
                            {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/Hsalinarum-pfeiffer2019-mobileElements.gff3",
                                type: "annotation",
                                searchable: true,
                                format: "gff3",
                                name: "IS annotation",
                                color: "#59A14F",
                                displayMode: "COLLAPSED"
                            },
                            
                            {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/mergedBRs_interaction-regions-entire-genome-rev.gff3",
                                type: "annotation",
                                format: "gff3",
                                name: "SmAP1 interaction (-)",
                                height: 30,
                                color: "#E15759",
                                displayMode: "SQUISHED"
                            },
                            
                            {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/TPS_gff_rev.gff3",
                                type: "annotation",
                                format: "gff3",
                                name: "TPS (-)",
                                height: 30,
                                color: "#59A14F",
                                displayMode: "SQUISHED"
                            },
                            
                            {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/ribosomal_RNA_TP2-rev.bedgraph.gz",
                                type: "wig",
                                format: "bedgraph",
                                autoscaleGroup: "ribo",
                                name: "Ribo-Seq TP2 (-)",
                                color: "darkgrey"
                            },
                            {
                                url: "https://alanlorenzetti.github.io/miscfiles/igv_tlr/total-RNA-TP2-fwd.bedgraph.gz",
                                type: "wig",
                                format: "bedgraph",
                                autoscaleGroup: "tot",
                                name: "RNA-Seq TP2 (-)",
                                color: "darkgrey"
                            }
                       ]
                    };
        
                  var igvDiv = document.getElementById("igvDiv");
        
                  igv.createBrowser(igvDiv, options)
                       .then(function (browser) {
                       console.log("Created IGV browser");
                       window.igvBrowser = browser;
                       })
                           
                  </script>
                  '
                  )
                )
            )
        ),
        tabItem(
            tabName = "help",
            fluidRow(
                box(
                    title = "Help",
                    width = 6, solidHeader = TRUE, status = "primary",
                    div(class = "my-class",
                        HTML("<h2>About the <i>Halobacterium salinarum</i> NRC-1 atlas</h2>"),
                        HTML("<p>
                        This resource is intended to help scientists on the manual inspection of gene sets displaying features 
                        that could explain underlying post-transcriptional regulation mechanisms. We developed this ShinyApp using
                        <a href=https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html>ComplexHeatmap</a> and
                        <a href=http://www.bioconductor.org/packages/release/bioc/html/InteractiveComplexHeatmap.html>InteractiveComplexHeatmap</a>
                        packages.
                        </p>"),
                        
                        HTML("<h2>How to use the atlas?</h2>"),
                        
                        HTML("<h3>Basic usage</h3>"),
                        HTML("<p>
                        There are two major ways to use the atlas:
                        <ol>
                        <li>Drag and click over the original heatmap to select regions of interest.
                        The selected region will become highlighted in the sub-heatmap panel;
                        </li>
                        <li>Click on the magnifier button at the original heatmap's bottom and type the gene names 
                        of interest. After clicking the \"Search\" button, the gene will become highlighted in the
                        sub-heatmap panel. One can also subset multiple genes by using a comma-separated list of genes.
                        For example: <code>VNG_1332G,VNG_7024,VNG_7025,VNG_7026,VNG_7039,VNG_7103</code>
                        </li>
                        </ol>
                        </p>
                        <p>The user can resize both panels by dragging the bottom-right corners.
                        Additionally, one can set up custom dimensions by clicking on the resize button
                        at the bottom of the panels and filling the fields.
                        </p>"),
                        
                        HTML("<h3>Intermediate usage</h3>"),
                        HTML("<p>
                        The user may perform gene selection using regular expressions as well. For that,
                        click on the magnifier button and tick the regular expression box. Type:
                        </p>
                        <p>
                        <code>VNG_702[5-8]</code>
                        </p>
                        <p>
                        This is equivalent to typing <code>VNG_7025,VNG_7026,VNG_7027,VNG_7028</code>
                        and should bring up the genes <i>gvpA</i>, <i>gvpC</i>, <i>gvpN</i>, and <i>gvpO</i>.
                        Most of genes have a code letter as the last character (e.g., VNG_1496G). In order to select those
                        genes using regular expression, it is a requirement to type the wildcard <code>.</code> at the end.
                        For example:
                        </p>
                        <p>
                        <code>VNG_149[0-9].</code>
                        </p>
                        <p>
                        Additional information about the available genes can be retrieved from our
                        <a href=https://alanlorenzetti.github.io/halo_nr_tx/>
                        <i>Halobacterium salinarum</i> non-redundant transcriptome resource.
                        </a>
                        </p>
                        <p>
                        It is possible to convert the sub-heatmap into a table
                        to observe the raw values. For that, find the table button (second tab) at
                        the bottom of the sub-heatmap. Then, click on the \"Open table\" button. It 
                        is even possible to customize the table by selecting the desired features in
                        the configure tab of the sub-heatmap (first tab).
                        </p>"),
                        
                        HTML("<h2>Static version</h2>"),
                        HTML("<p>
                        For those who are not familiar with interactive heatmaps or those who
                        are experiencing technical difficulties, we made available a static
                        version of this atlas in PDF. Click <a href=\"data/abundanceHeatmap_expanded_en.pdf\">here</a>
                        to open the file.
                        </p>")
                    )
                )
            )
        )
    )
)

# Define UI for application that draws a histogram
ui <- dashboardPage(
    header = header,
    sidebar = sidebar,
    body = body)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    makeInteractiveComplexHeatmap(input, output, session, htComplete, "ht_atlas")
    
    output$heatmapLegend = renderImage({
        list(src = "data/heatmap_legends.png",
             width = 1000,
             height = 500)
    }, deleteFile = F)
    
    # #getting igv shiny to work
    # output$igvShiny_0 <- renderIgvShiny(
    #     igvShiny(list(
    #         genomeName="hsal",
    #         initialLocus="VNG_1496G",
    #         displayMode="SQUISHED")
    #     )
    # )
    
    addResourcePath("data", "data")
}

# Run the application 
shinyApp(ui = ui, server = server)

