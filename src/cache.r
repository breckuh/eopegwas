cache.init = function (path = "cache/", verbose = TRUE) {
  if (verbose & exists("cache._storage")) {
    print("Notice: cache already initialized. Keeping existing cache.")
  } else {
    cache._path <<- path
    cache._storage <<- list()  
  }
}

cache.clearAll = function () {
  lapply(names(cache._storage), cache.clear)
}

cache.clear = function(key) {
  if (!cache.has(key)) {
    warning(glue("Cache key {key} does not exist. Not clearing anything."))
  } else {
    cache._storage[[key]] <<- NULL
    file.remove(cache._makeLink(key))
  }
}

cache._makeLink = function (key) {
  glue("{cache._path}{key}.rds")
}

cache.has = function (key) {
  if (!is.null(cache._storage[[key]]))
    return(TRUE)
  outFile = cache._makeLink(key)
  if (file.exists(outFile))
    return(TRUE)
  return(FALSE)
}

cache.get = function (key, dataFn, reset = FALSE) {
  if (!is.null(cache._storage[[key]]))
    return(cache._storage[[key]])
  outFile = cache._makeLink(key)
  if (reset | !file.exists(outFile)) {
    value = dataFn() 
    saveRDS(value, outFile)
  } else {
    value = readRDS(outFile)
  }
  cache._storage[[key]] <<- value
  return(value)
}