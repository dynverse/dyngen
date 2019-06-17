#include <Rcpp.h>
#include "ssa.hpp"
#include "ssa_em.hpp"
#include "ssa_direct.hpp"
#include "utils.hpp"

using namespace Rcpp;

typedef void (*TR_FUN)(const NumericVector&, const NumericVector&, const NumericVector&, NumericVector&);

// [[Rcpp::export]]
SEXP make_ssa_em(double h, double noise_strength) {
  SSA_EM *ssa = new SSA_EM(h, noise_strength);
  XPtr<SSA_EM> ptr(ssa);
  return ptr;
}

// [[Rcpp::export]]
SEXP make_ssa_direct() {
  SSA_direct *ssa = new SSA_direct();
  XPtr<SSA_direct> ptr(ssa);
  return ptr;
}

// [[Rcpp::export]]
List simulate(
    SEXP transition_fun,
    SEXP ssa_alg,
    const NumericVector& initial_state,
    const NumericVector& params,
    const NumericMatrix& nu,
    const double final_time,
    const double max_duration,
    const bool stop_on_neg_state,
    const bool stop_on_neg_propensity,
    const bool verbose,
    const double console_interval
) {
  List output(10);
  
  TR_FUN transition_fun_ = *XPtr<TR_FUN>(transition_fun);
  SSA *ssa_alg_ = XPtr<SSA>(ssa_alg);
  
  NumericVector state = clone(initial_state);
  NumericVector sim_time = {0};
  NumericVector transition_rates(nu.ncol());
  
  // calculate initial transition rates
  transition_fun_(state, params, sim_time, transition_rates);
  
  output(0) = List::create(
    _["time"] = sim_time,
    _["state"] = clone(state),
    _["transition_rates"] = clone(transition_rates)
  );
  int output_nexti = 1;
  
  int realtime_start = time(NULL);
  int realtime_nextconsole = realtime_start, realtime_nextinterrupt = realtime_start, realtime_curr = realtime_start;
  
  if (verbose) {
    Rcout << "Running SSA with console output every " << console_interval << " seconds" << std::endl;
    Rcout << "Start time: " << "CURRTIME" << std::endl;
  }
  /*if (verbose) {
   cat("Running ", method$name, " method with console output every ", console.interval, " time step\n", sep = "")
   cat("Start wall time: ", format(time.start), "\n" , sep = "")
   utils::flush.console()
  }*/
  
  while (sim_time[0] < final_time && (realtime_curr - realtime_start) <= max_duration)  {
    realtime_curr = time(NULL);
    
    if (realtime_nextinterrupt <= realtime_curr) {
      checkUserInterrupt();
      realtime_nextinterrupt += 1;
    }
    
    if (verbose && realtime_nextconsole <= realtime_curr) {
      Rcout << realtime_curr << " | sim_time = " << sim_time[0] << " : " << "STATE" << std::endl;
      realtime_nextconsole += console_interval;
    }
    
    NumericVector dtime(1);
    NumericVector dstate(state.size());
    // TODO: pass time and state directly to step, instead of dtime and dstate?
    
    ssa_alg_->step(state, transition_rates, nu, dtime, dstate);
    // step(state, transition_rates, nu, dtime, dstate);
    
    state += dstate;
    sim_time[0] += dtime[0];
    
    // Rcout << "sim_time: " << sim_time << ", dtime: " << dtime << std::endl;
    
    
    /*# Check that no states are negative (can occur in some tau-leaping methods)
     invalid_ix <- is.na(state) | state < 0
     if (any(invalid_ix)) {
     if (stop_on_neg_state) {
     stop("state vector contains negative values\n", paste(names(state)[invalid_ix], collapse = ", "))
     } else {
     state[invalid_ix] <- 0
     }
     }*/
    
    transition_fun_(state, params, sim_time, transition_rates);
    
    /*
     invalid_ix <- is.na(transition_rates) | transition_rates < 0
     if (any(invalid_ix)) {
     if (stop_on_neg_propensity) {
     stop("transition rate contains negative values\n", paste(names(transition_rates)[invalid_ix], collapse = ", "))
     } else {
     transition_rates[invalid_ix] <- 0
     }
     }
     */
    
    if (output_nexti == output.size()) {
      output = resize(output, output.size() * 2);
    }
    
    output(output_nexti) = List::create(
      _["time"] = sim_time,
      _["state"] = clone(state),
      _["transition_rates"] = clone(transition_rates)
    );
    output_nexti++;
  }
  
  // Display the last time step on the console
  if (verbose) {
    Rcout << "sim_time = " << sim_time[0] << " : " << "STATE" << std::endl;
  }
  
  /*
   int realtime_end = time(NULL);
   int realtime_elapsed = realtime_end - realtime_start;
   
   DataFrame stats = DataFrame::create(
   // _["method"] = name,
   _["final_time_reached"] = sim_time >= final_time,
   // _["extinction"] = all(state == 0),
   // _["negative_state"] = any(state < 0),
   // _["zero_prop"] = all(transition_rates == 0),
   _["max_duration"] = realtime_elapsed >= max_duration,
   _["realtime_start"] = realtime_start,
   _["realtime_end"] = realtime_end,
   _["realtime_elapsed"] = realtime_elapsed,
   );
   // if (verbose) {
   //   //print(stats)
   // }
   
   return List::create(
   _["output"] = output,
   _["stats"] = stats
   );
   */
  return(output);
}
