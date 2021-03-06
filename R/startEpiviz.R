.register_all_the_epiviz_things <- function(app) {
  # register actions requested from epiviz app
  app$server$register_action("getMeasurements", function(request_data) {
    app$data_mgr$get_measurements()
  })
  
  app$server$register_action("getRows", function(request_data) {
    app$data_mgr$get_rows(request_data$seqName,
      request_data$start,
      request_data$end,
      request_data$metadata, 
      request_data$datasource)
  })
  
  app$server$register_action("getValues", function(request_data) {
    app$data_mgr$get_values(request_data$seqName, 
      request_data$start, 
      request_data$end,
      request_data$datasource, 
      request_data$measurement)
  })
  
  app$server$register_action("getSeqInfos", function(request_data) {
    app$data_mgr$get_seqinfo()
  })
    
  app$server$register_action("registerChartTypes", function(request_data) {
    app$chart_mgr$.register_available_chart_types(request_data$data)
  })
  
  app$server$register_action("setChartSettings", function(request_data) {
    app$chart_mgr$.update_chart_settings(request_data)
  })
  
  ## TODO: register action 'search'
}

#' Start epiviz app and create \code{\link{EpivizApp}} object to manage connection.
#' 
#' @param host (character) use a host for the epiviz app other than the cbcb.umd.edu hosts.
#' @param http_port (integer) port at host serving the epiviz app.
#' @param path (character) path at host where epiviz app is located.
#' @param use_devel (logical) use the devel epiviz application server (http://epiviz-dev.cbcb.umd.edu).
#' @param chr (character) chromosome to browse to on app startup.
#' @param start (integer) start location to browse to on app startup.
#' @param end (integer) end location to browse to on app startup.
#' @param debug (logical) start epiviz app in debug mode.
#' @param workspace (character) a workspace id to load in the epiviz app on startup.
#' @param scripts (character) URLs for JavaScript plugin scripts to be imported when epiviz is loaded (see \url{http://epiviz.cbcb.umd.edu/help} for details).
#' @param gists (character) Ids for github gists (\url{http://gist.github.com}) containing JavaScript plugin scripts to
#'  be imported when epiviz is loaded (see \url{http://epiviz.cbcb.umd.edu/help} for details).
#' @param use_cookie (logical) use cookies within the epiviz app.
#' @param register_function (function) function used to register actions and charts on the epiviz app.
#' @param open_browser (logical) browse to the epiviz URL before exiting function.
#' @param server (EpivizServer) if not \code{NULL} use this object as underlying WebSocket and HTTP server
#' @param browser_fun (function) function used to browse to URL (\code{browseURL} by default)
#' @param ws_host (character) host address to use for websocket connection ("localhost" by default)
#' @param ... additional parameters passed to \code{\link[epivizrServer]{createServer}}.
#' 
#' @return An object of class \code{\link{EpivizApp}}
#' 
#' @examples
#' # see package vignete for example usage
#' app <- startEpiviz(non_interactive=TRUE, open_browser=TRUE)
#' app$stop_app()
#' 
#' @export
startEpiviz <- function(host=NULL, http_port=NULL, path=NULL, use_devel=FALSE, 
                        chr="chr11", start=99800000, end=103383180, 
                        debug=FALSE, workspace=NULL, scripts=NULL, gists=NULL, use_cookie=TRUE,
                        register_function=.register_all_the_epiviz_things,
                        open_browser=TRUE, server=NULL, browser_fun=utils::browseURL, 
                        ws_host="localhost", ...) {

  if(is.null(server) || !is(server, "EpivizServer") ) {
    server <- epivizrServer::createServer(...)
  }
  
  data_mgr <- epivizrData::createMgr(server)
  chart_mgr <- EpivizChartMgr$new(server)
  
  url_parms <- list(host=host,
    http_port=http_port,
    path=path,
    ws_port=server$.port,
    use_devel=use_devel,
    debug=debug,
    chr=chr,
    start=start,
    end=end,
    workspace=workspace,
    scripts=scripts,
    gists=gists,
    use_cookie=use_cookie,
    ws_host=ws_host)
  
  app <- EpivizApp$new(.url_parms=url_parms,
                       .browser_fun=browser_fun,
                       server=server,
                       data_mgr=data_mgr,
                       chart_mgr=chart_mgr)
  
  if (app$server$.verbose) {
    cat("Starting Epivizr!\n")
  }
  
  register_function(app)
  
  if (!app$server$is_interactive()) {
    return(app)
  }
  
  if (app$server$.verbose) {
    cat("Initializing session manager...\n")
  }
  tryCatch({
    if (app$server$.verbose) {
      cat("Opening connections...\n")
    }
    
    if (open_browser) {
      app$.open_browser()
    }

    if (app$server$.verbose) {
      cat("Done starting Epivizr!\n")
    }
    app
  }, error=function(e) {
    app$stop_app()
    stop("Error starting Epiviz: ", e)
  }, interrupt=function(e) {NULL})  
}
