##
## ---- Record the installed package build used by a throughput benchmark:
##
fn_benchmark_package_build_metadata <-  function(package_name = "BayesMVP") {
      package_directory <-  find.package(package_name)
      package_libraries <-  list.files(file.path(package_directory, "libs"),
                                      pattern = paste0("\\", .Platform$dynlib.ext, "$"),
                                      recursive = TRUE,
                                      full.names = TRUE)
      list(package = package_name,
           version = as.character(utils::packageVersion(package_name)),
           library_path = package_directory,
           library_md5 = tools::md5sum(package_libraries),
           description_md5 = tools::md5sum(file.path(package_directory, "DESCRIPTION")),
           R_code_md5 = tools::md5sum(list.files(file.path(package_directory, "R"), full.names = TRUE)),
           R_version = R.version.string)
}

##
## ---- Save timings with their settings; refuse to replace a different experiment:
##
fn_save_benchmark_timing_results <-  function(timing_array,
                                             file_path,
                                             benchmark_settings) {
      if (file.exists(file_path)) {
            previous_settings <-  attr(readRDS(file_path), "benchmark_settings")
            if (!identical(previous_settings, benchmark_settings)) {
                  stop("Existing benchmark has different or unrecorded settings: ", file_path,
                       ". Choose a new benchmark_results_subdirectory to preserve both experiments.")
            }
      }
      dir.create(dirname(file_path), recursive = TRUE, showWarnings = FALSE)
      attr(timing_array, "benchmark_settings") <-  benchmark_settings
      attr(timing_array, "saved_at") <-  Sys.time()
      saveRDS(object = timing_array, file = file_path)
      invisible(file_path)
}






















