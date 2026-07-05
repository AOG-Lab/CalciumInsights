# File import helper used by both active modules.
load_file <- function(name, path, ext) {
  ext <- tools::file_ext(name)
  switch(ext,
         csv = vroom::vroom(path, delim = ","),
         tsv = vroom::vroom(path, delim = "\t"),
         json = jsonlite::fromJSON(path),  # Suponiendo que "jisonreader" es la biblioteca adecuada
         stop("Formato de archivo no compatible.")
  )
}
