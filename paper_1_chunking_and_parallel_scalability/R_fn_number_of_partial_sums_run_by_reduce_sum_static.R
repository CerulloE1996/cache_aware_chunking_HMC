##
## ======================================================================================================================================
## R_fn_number_of_partial_sums_run_by_reduce_sum_static.R
##
## The number of partial sums which Stan's reduce_sum_static() actually runs for a requested N_chunks.
##
## The Stan models of Paper 1 (AD_Stan_tape_chunked, AD_Stan_WCP and AD_Stan_WCP_chunking; the model file
## LC_MVP_bin_PartialLog_v5_reduce_sum_static.stan) pass chunk_size = ceiling(N / N_chunks) to reduce_sum_static() as its
## grainsize. reduce_sum_static() runs tbb::parallel_reduce() with tbb::simple_partitioner on a tbb::blocked_range of the N
## individuals (stan/math/prim/functor/reduce_sum_static.hpp), which splits a range at its midpoint (the first part holds
## floor(size / 2) individuals) for as long as it holds more than grainsize individuals. The number of partial sums is the
## number of leaves of this recursive halving. For every N_chunks tested in Paper 1 this is the next power of two
## (e.g. N_chunks = 10 runs as 16 partial sums at N = 500, and N_chunks = 250 as 256 at N = 50,000), but it is not a power of
## two in general (e.g. N = 5 with N_chunks = 3 runs as 3 partial sums).
##
## The container-chunked Stan model (AD_Stan_chunked) loops over ceiling(N / chunk_size) blocks itself, and NicoStan+BayesMVP
## chunks in its own loop, so neither is affected.
##
#### ---- the number of leaves of the recursive halving of N_units with grainsize ceiling(N_units / N_chunks_requested) ----------------
#' The number of partial sums which reduce_sum_static() runs for a requested N_chunks
#'
#' @param N_units The number of units which reduce_sum_static() slices (the N individuals of the LC-MVP model).
#' @param N_chunks_requested The requested N_chunks, i.e. grainsize = ceiling(N_units / N_chunks_requested). NA stays NA.
#' @return The number of partial sums, one for each element of N_units and N_chunks_requested (recycled).
fn_number_of_partial_sums_run_by_reduce_sum_static <-  function( N_units,
                                                                 N_chunks_requested
) {

        fn_number_of_partial_sums_for_one_case <-  function( N_units_one,
                                                             N_chunks_requested_one
        ) {

                if (is.na(x = N_units_one) || is.na(x = N_chunks_requested_one)) return(NA_real_)
                ##
                grainsize <-  ceiling(N_units_one / N_chunks_requested_one)
                ##
                fn_number_of_leaves <-  function(N_units_in_range) {

                        if (N_units_in_range <= grainsize) return(1)
                        ##
                        N_units_in_first_part <-  N_units_in_range %/% 2
                        ##
                        return(fn_number_of_leaves(N_units_in_first_part) + fn_number_of_leaves(N_units_in_range - N_units_in_first_part))

                }
                ##
                return(fn_number_of_leaves(N_units_one))

        }
        ##
        return(unname(mapply(FUN = fn_number_of_partial_sums_for_one_case, N_units, N_chunks_requested)))

}
























