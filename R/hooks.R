## KNITR HOOK

epiviz.chart <- function(before, options, envir) {
  if (!before) {
    charts <- options$chart.names
    if (is.null(charts)) {
      stop(options$label, " must have an option, chart.names, with list of variable names of charts")
    }
    
    chart_accumulator <- vector("list", length(charts))
    
    for (i in 1:length(charts)) {
      chart_obj <- envir[[charts[[i]]]]
      if(is.null(chart_obj)) {
        stop(charts[[i]], " doesn't exist in environment")
      }

      row_data <- row_data(chart_obj)
      col_data <- NULL
      # Blocks Tracks and Genes Tracks do not use values
      if (chart_obj$.type != "epiviz.plugins.charts.BlocksTrack" && 
          chart_obj$.type != "epiviz.plugins.charts.GenesTrack") {
        col_data <- col_data(chart_obj)                  
      }
      
      data <- paste0(row_data, col_data, sep=",")
      poly_chart <- chart_obj$to_polymer(data)
      
      chart_accumulator[[i]] <- poly_chart  
    }
    
    # If chunk option was set with epivizChart=list(epivizEnvir=TRUE,...)
    # if (!is.null(options$epiviz.envir) && options$epiviz.envir) {
    #   charts <- paste0(chart_accumulator, collapse='\n')
    #   
    #   epivizEnvir <- paste0(epiviz_environment(before=TRUE, app$.url_parms$chr,
    #     app$.url_parms$start, app$.url_parms$end), charts, epiviz_environment(before=FALSE))
    #   
    #   return(epiviznEnvir)
    # } else {
      # Chart(s) will not be wrapped in an epiviz environment
      epivizChunkCharts <- paste0(chart_accumulator, collapse='\n')
    }
    
    return(epivizChunkCharts)
}


row_data <- function(chart_obj) {
  ms_obj <- app$get_ms_object(chart_obj)
  # TODO: Fix calling  global variable 'app'
  query <- GenomicRanges::GRanges(app$.url_parms$chr, ranges = IRanges::IRanges(app$.url_parms$start, app$.url_parms$end))
  
  result <- ms_obj$get_rows(query = query)
  json_row_data <- epivizrServer::json_writer(result)
  
  return(json_row_data)
}

col_data <- function(chart_obj) {
  ms_obj <- app$get_ms_object(chart_obj)
  # TODO: Fix calling  global variable 'app'
  query <- GenomicRanges::GRanges(app$.url_parms$chr, ranges = IRanges::IRanges(app$.url_parms$start, app$.url_parms$end))
  
  ms_list <- ms_obj$get_measurements()
  col_data <- vector("list", length(ms_list)) 
  
  for(i in 1:length(ms_list)) {
    ms <- ms_list[[i]]
    values <- ms_obj$get_values(query=query, measurement=ms@id)
    values_json <- epivizrServer::json_writer(values)
    
    col_data[[i]] <- values_json
  }
  
  col_data <- paste0(col_data, collapse='')
  return(col_data)
}

# Data source
# epivizDS <- function(before, options, envir) {
#   if (!before) {
#     for (obj in ls(envir)) {
#       app <- eval(parse(text=obj))
#       
#       if (is(app, "EpivizApp")) {
#         provider_url <- paste0("ws://", app$.url_parms$ws_host, ":", app$.url_parms$ws_port)
# 
#         epivizDataSource <- epiviz_data_source(
#           provider_type="epiviz.data.WebsocketDataProvider",
#           provider_id="app", 
#           provider_url=provider_url
#         )
#         return(epivizDataSource)
#       }
#     }
#   }
# }


####################
# POLYMER ##########
####################

# epiviz_data_source <- function(provider_type, provider_id, provider_url){
#   epiviz_data_source <- paste0(
#     '<epiviz-data-source ',
#     'provider-type="', provider_type, '"',
#     'provider-id="', provider_id, '"',
#     'provider-url="', provider_url, '"',
#     '>',
#     '</epiviz-data-source>'
#   )
#   
#   return(epiviz_data_source)
# }
# 
# epiviz_environment <- function(before=TRUE, chr="chr11", start=99800000, end=103383180) {
#   epiviz_environment <- NULL 
#   
#   if (before) {
#     epiviz_environment <- paste0(
#       '<epiviz-environment ', 
#       'chr="', chr,'\"',
#       'start="', start, '"',
#       'end="', end, '"',
#       '>'
#     )
#   } else {
#     epiviz_environment <- "</epiviz-environment>"
#   }
#   
#   return(epiviz_environment)
# }
# 
# polymer_chart_type <- function(chart_type) {
#   if (chart_type == "epiviz.plugins.charts.BlocksTrack") return("epiviz-blocks-track")
#   if (chart_type == "epiviz.plugins.charts.GenesTrack")  return("epiviz-genes-track")
#   if (chart_type == "epiviz.plugins.charts.HeatmapPlot") return("epiviz-heatmap-plot")
#   if (chart_type == "epiviz.plugins.charts.LinePlot") return("epiviz-line-plot")
#   if (chart_type == "epiviz.plugins.charts.LineTrack") return("epiviz-line-track")
#   if (chart_type == "epiviz.plugins.charts.ScatterPlot") return("epiviz-scatter-plot")
#   if (chart_type == "epiviz.plugins.charts.StackedLinePlot") return("epiviz-stacked-line-plot")
#   if (chart_type == "epiviz.plugins.charts.StackedLineTrack") return("epiviz-stacked-line-track")
# }
# polymer_chart_type_json <- function(chart_type) {
#   if (chart_type == "epiviz.plugins.charts.BlocksTrack") return("epiviz-json-blocks-track")
#   if (chart_type == "epiviz.plugins.charts.GenesTrack")  return("epiviz-json-genes-track")
#   if (chart_type == "epiviz.plugins.charts.HeatmapPlot") return("epiviz-json-heatmap-plot")
#   if (chart_type == "epiviz.plugins.charts.LinePlot") return("epiviz-json-line-plot")
#   if (chart_type == "epiviz.plugins.charts.LineTrack") return("epiviz-json-line-track")
#   if (chart_type == "epiviz.plugins.charts.ScatterPlot") return("epiviz-json-scatter-plot")
#   if (chart_type == "epiviz.plugins.charts.StackedLinePlot") return("epiviz-json-stacked-line-plot")
#   if (chart_type == "epiviz.plugins.charts.StackedLineTrack") return("epiviz-json-stacked-line-track")
# }
# 
# polymer_chart <- function(chart_type, chart_id, measurements, data) {
#   polymer_chart <- paste0(
#     '<', chart_type, ' ',
#     'class="charts" ',
#     'id="', chart_id, '"',
#     "measurements='", measurements, "'",
#     "data='", data, "'",
#     '>',
#     '</', chart_type, '>', "\n"
#   )
#   
#   return(polymer_chart)
# }
# 




knitr::knit_hooks$set(
  # epivizDS=epivizDS, 
  epiviz.chart=epiviz.chart)


