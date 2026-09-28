ws_dir <-  getwd()   ## run from working_set/stan
so_path <-  "LC_MVP_bin_PartialLog_v5_reduce_sum_static_model.so"   ## the model compiled with BridgeStan, as in the benchmarks (see build_and_run.sh)
for (case in c("N500_cs500", "N20000_cs5000")) {
    N <-  sub(pattern = "_cs.*", replacement = "", x = sub(pattern = "^N", replacement = "", x = case))
    model <-  bridgestan::StanModel$new(lib = so_path, data = file.path(ws_dir, "json", paste0("data_", case, ".json")), seed = 0)
    for (pt in c("A", "B")) {
        init_json <-  paste(readLines(con = file.path(ws_dir, "json", paste0("init_", pt, "_N", N, ".json"))), collapse = "\n")
        theta <-  model$param_unconstrain_json(init_json)
        out <-  model$log_density_gradient(theta_unc = theta, propto = TRUE, jacobian = TRUE)
        message(paste0("BRIDGESTAN ", pt, "_", case, " P=", length(theta), " lp=", formatC(out$val, digits = 10, format = "g"),
                       " grad_l2=", formatC(sqrt(sum(out$gradient^2)), digits = 10, format = "g")))
    }
}
