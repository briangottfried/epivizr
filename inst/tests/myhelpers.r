setEpivizrTestOpts <- function(sendRequest=TRUE,
                                daemonized=TRUE,
                               local=FALSE,
                                devel=FALSE,
                               test=FALSE,
                                debug=TRUE,
                               port=7312L) {
  url <- if (devel) "epiviz-dev" else "epiviz"

  if (!local) {
    url <- sprintf("http://%s.cbcb.umd.edu", url)
  } else {
    url <- sprintf("http://localhost/~hcorrada/%s", url)
  }
  
  if (test) {
    url <- paste(url, "test_socket.php", sep="/")
  } else {
    url <- paste(url, "index.php", sep="/")
  }

  options(epivizrTestSendRequest=sendRequest,
          epivizrTestDaemonized=daemonized,
          epivizrTestURL=url,
          epivizrTestDebug=debug,
          epivizrTestPort=port)
  invisible()
}

getEpivizrTestOpts=function() {
  out <- list(sendRequest=getOption("epivizrTestSendRequest"),
              daemonized=getOption("epivizrTestDaemonized"),
              url=getOption("epivizrTestURL"),
              debug=getOption("epivizrTestDebug"),
              port=getOption("epivizrTestPort"))
  print(out)
}

setEpivizrTestOpts()

test_srv=function(dem=TRUE) {setEpivizrTestOpts(daemonized=dem); test(filter=".*server.*")}
test_reg=function() test(filter=".*register.*")
test_mes=function(req=TRUE,dem=TRUE) {setEpivizrTestOpts(sendRequest=req, daemonized=dem); test(filter=".*Measure.*")}
test_fet=function(req=TRUE) {setEpivizrTestOpts(sendRequest=req); test(filter=".*fetch.*")}
test_cha=function(req=TRUE,dem=TRUE) {setEpivizrTestOpts(sendRequest=req, daemonized=dem); test(filter=".*Charts.*")}
test_dev=function(req=TRUE,dem=TRUE) {setEpivizrTestOpts(sendRequest=req, daemonized=dem); test(filter=".*Device.*")}

test_some=function(req=TRUE,dem=TRUE) {test_mes(req=req,dem=dem); test_cha(req=req,dem=dem); test_dev(req=req,dem=dem)}

testb=function() {test_srv(FALSE);test_reg();test_fet(FALSE)}
testb1=function() {test_srv(TRUE); test_reg(); test_fet(TRUE)}
test0=function() test_some(FALSE,FALSE)
test1=function() test_some(TRUE,FALSE)
test2=function() test_some(TRUE,TRUE)

test_all=function() {testb(); testb1(); test0(); test1(); test2()}

