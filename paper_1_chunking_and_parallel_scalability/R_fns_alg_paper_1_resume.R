##
## ======================================================================================================================================
## R_fns_alg_paper_1_resume.R
##
## Cache helpers for the unified Paper 1 study. A saved result is reused when the case it measured is the same case:
## device, algorithm, backend, N, chunks, threads, chains, threads per chain, iterations, timing method, repeat, seeds and
## Mplus iteration mode. Nothing else (model file contents, package builds, runtime records, settings fingerprints) is
## compared. Superseded duplicate measurements are kept outside the cases folder, so each case has one saved result.
##

fn_paper1_resume_object_md5 <-  function(object) {

        temporary_file <-  tempfile(pattern = "paper1_resume_signature_")
        on.exit(unlink(x = temporary_file), add = TRUE)
        saveRDS(object = object, file = temporary_file, compress = FALSE, version = 2)
        unname(tools::md5sum(files = temporary_file))

}
##
fn_paper1_resume_scalar <-  function(value) {

        if (is.null(value)) return(NULL)
        if (length(value) == 1L && is.atomic(value)) {
            if (is.numeric(value)) return(as.numeric(value))
            if (is.integer(value)) return(as.numeric(value))
            if (is.logical(value)) return(as.logical(value))
            if (is.character(value)) return(as.character(value))
        }
        value

}
##
fn_paper1_resume_named_list <-  function(value) {

        if (is.null(value)) return(NULL)
        if (!is.list(value)) return(fn_paper1_resume_scalar(value))
        value <- lapply(value, fn_paper1_resume_named_list)
        if (!is.null(names(value))) value <- value[sort(names(value), method = "radix")]
        value

}
##
fn_paper1_resume_case_work <-  function(case) {

        required_fields <- c("device", "algorithm", "N", "num_chunks", "n_threads", "n_chains",
                             "threads_per_chain", "n_iter", "run", "seed")
        missing_fields <- setdiff(required_fields, names(case))
        if (length(missing_fields)) {
            stop("Cannot form a Paper 1 resume key; case is missing: ", paste(missing_fields, collapse = ", "))
        }
        ## These fields alter the sampler calls or the short/long raw pair.
        ## case_id, status, error, timing_estimator, and all measured outputs
        ## are intentionally excluded.  timing_estimator is a reporting choice
        ## and may be changed while the saved raw pair remains reusable.
        work_fields <- c("device", "algorithm", "execution_backend", "N", "num_chunks", "n_threads",
                         "n_chains", "threads_per_chain", "n_iter", "n_iter_short_run", "timing_method",
                         "run", "base_seed", "seed",
                         "mplus_iteration_mode", "mplus_iteration_mode_short_run", "stan_chunk_size")
        work <- lapply(work_fields, function(field) {
            value <- case[[field]]
            if (is.null(value)) return(NULL)
            if (length(value) != 1L) return(value)
            fn_paper1_resume_scalar(value)
        })
        names(work) <- work_fields
        ## Historical/current grid rows can omit this derived Stan field.  It is
        ## part of the actual reduce_sum partition and must therefore be keyed.
        if (is.null(work$stan_chunk_size) &&
            as.character(case$algorithm) %in% c("AD_Stan_chunked", "AD_Stan_tape_chunked", "AD_Stan_WCP")) {
            work$stan_chunk_size <- as.numeric(ceiling(as.numeric(case$N) / as.numeric(case$num_chunks)))
        }
        work

}
##
fn_paper1_resume_key <-  function(case, settings = NULL, dataset_md5 = NULL, metadata = NULL) {

        fn_paper1_resume_object_md5(fn_paper1_resume_named_list(list(schema_version = 2L,
                                                                     case = fn_paper1_resume_case_work(case = case))))

}
##
fn_paper1_resume_cache_path <-  function(cache_dir, key) {

        if (length(key) != 1L || is.na(key) || !grepl("^[[:xdigit:]]{32}$", key)) return(NA_character_)
        file.path(cache_dir, paste0("case_", key, ".rds"))

}
##
fn_paper1_resume_validate_row <-  function(row, ...) {

        is.data.frame(row) && nrow(row) == 1L &&
            "status" %in% names(row) && identical(as.character(row$status[[1]]), "completed") &&
            "elapsed_seconds" %in% names(row) && isTRUE(is.finite(as.numeric(row$elapsed_seconds[[1]])))

}
##
fn_paper1_resume_cache_index <-  function(cache_dir, settings = NULL) {

        index <- new.env(parent = emptyenv())
        files <- sort(list.files(cache_dir, pattern = "^case_[[:xdigit:]]{32}\\.rds$", full.names = TRUE), method = "radix")
        for (file in files) {
            stored <- tryCatch(readRDS(file), error = function(error) NULL)
            if (!is.list(stored) || !isTRUE(fn_paper1_resume_validate_row(stored$row))) next
            key <- tryCatch(fn_paper1_resume_key(case = stored$row), error = function(error) NA_character_)
            if (length(key) != 1L || is.na(key)) next
            if (!exists(key, envir = index, inherits = FALSE)) assign(key, stored, envir = index)
        }
        index

}
##
fn_paper1_cache_read <-  function(cache_dir, key, index = NULL, ...) {

        if (!is.null(index) && exists(key, envir = index, inherits = FALSE)) return(get(key, envir = index)$row)
        cache_file <- fn_paper1_resume_cache_path(cache_dir = cache_dir, key = key)
        if (is.na(cache_file) || !file.exists(cache_file)) return(NULL)
        stored <- tryCatch(readRDS(cache_file), error = function(error) NULL)
        if (!is.list(stored) || !isTRUE(fn_paper1_resume_validate_row(stored$row))) return(NULL)
        stored$row

}
##
fn_paper1_cache_write <-  function(cache_dir,
                                   key,
                                   row,
                                   metadata = NULL,
                                   preserve_historical_stan_runs = integer()
) {

        cache_file <- fn_paper1_resume_cache_path(cache_dir = cache_dir, key = key)
        if (is.na(cache_file)) return(invisible(FALSE))
        if (!fn_paper1_resume_validate_row(row = row,
                                           preserve_historical_stan_runs = preserve_historical_stan_runs)) stop("Only one completed Paper 1 result row can be cached.")
        dir.create(path = cache_dir, recursive = TRUE, showWarnings = FALSE)
        ## Preserve an existing valid completion for the same exact key.  This
        ## keeps reruns idempotent and avoids replacing a good row concurrently.
        if (file.exists(cache_file) && !is.null(fn_paper1_cache_read(cache_dir = cache_dir, key = key,
                                                                      preserve_historical_stan_runs = preserve_historical_stan_runs))) return(invisible(FALSE))
        temporary_file <- tempfile(pattern = paste0(".case_", key, "_"), tmpdir = cache_dir)
        on.exit(unlink(x = temporary_file), add = TRUE)
        stored <- list(schema_version = 1L, key = key, row = row, metadata = metadata,
                       stored_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
        saveRDS(object = stored, file = temporary_file, compress = FALSE, version = 2)
        if (!file.rename(from = temporary_file, to = cache_file)) {
            ## Some platforms refuse rename-over-existing.  Preserve the
            ## malformed entry for inspection, then install the replacement.
            if (file.exists(cache_file)) {
                invalid_file <- paste0(cache_file, ".invalid_", format(Sys.time(), "%Y%m%d_%H%M%S"))
                if (!file.rename(from = cache_file, to = invalid_file)) unlink(cache_file)
            }
            if (!file.rename(from = temporary_file, to = cache_file) && !file.exists(cache_file)) {
                stop("Could not atomically install Paper 1 cache file: ", cache_file)
            }
        }
        invisible(TRUE)

}






















