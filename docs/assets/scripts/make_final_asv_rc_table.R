make_final_asv_rc_table <- function(df, page_size = 5) {
  reactable(
    df, defaultColDef = colDef(
        header = function(value) gsub("_", " ", value, fixed = TRUE),
        cell = function(value) format(value, nsmall = 0),
        align = "center", 
        filterable = FALSE, 
        sortable = TRUE,
        resizable = TRUE, 
        footerStyle = list(fontWeight = "bold")
    ),
    columns = list(
        SampleID = colDef(
            name = "SampleID", 
            sticky = "left", 
            style = list(borderRight = "1px solid #eee"),
            headerStyle = list(borderRight = "1px solid #eee"),
            align = "left", 
            minWidth = 150
            ),
        raw_rc = colDef(
            name = "raw rc"
            ),
        cutadapt_rc = colDef(
            name = "rc after cutadapt"
            ),
        lotus3_rc = colDef(
            name = "rc after LotuS3"
            ),
        lotus3_asv = colDef(
            name = "ASVs after LotuS3"
            ),
        final_rc = colDef(
             name = "final rc"
             ),
        final_asv = colDef(
          name = "final ASVs"
             )
    ),
    searchable = FALSE, 
    defaultPageSize = page_size, 
    pageSizeOptions = c(5, 10, nrow(df)),
    showPageSizeOptions = TRUE, 
    highlight = TRUE, 
    bordered = TRUE,
    striped = TRUE, 
    compact = FALSE, 
    wrap = FALSE, 
    showSortable = TRUE,
    fullWidth = TRUE, 
    theme = reactableTheme(style = list(fontSize = "0.8em"))
    ) 
}