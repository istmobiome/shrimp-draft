# --- footer.R ---
## source_code("my_username", "my_repo")

source_code <- function(user_name, repo_name) {
  rm(list = ls(pattern = "tmp_"))
  
  tmp_in <- knitr::current_input(dir = FALSE)
  tmp_in <- xfun::with_ext(tmp_in, ".qmd")
  
  # Absolute path to current file
  tmp_abs <- normalizePath(tmp_in)
  
  # Instead of ".", use project root from Quarto
  proj_root <- rprojroot::find_root(rprojroot::has_file("_quarto.yml"))
  
  # Relative path from project root
  rel_path <- fs::path_rel(tmp_abs, start = proj_root)
  
  # Construct GitHub URL
  tmp_ghub <- paste0("https://github.com/", user_name, "/", repo_name, "/blob/main/", rel_path)
  
  return(tmp_ghub)
}

page_name <- function() {
  rm(list = ls(pattern = "tmp_"))
  tmp_in <- knitr::current_input(dir = FALSE)
  tmp_in <- xfun::with_ext(tmp_in, ".qmd")
  return(tmp_in)                                 
}

get_time <- function() {
  Sys.Date()
}

