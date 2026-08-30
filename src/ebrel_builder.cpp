//------------------------------ Ebrel builder ---------------------------------

#include "ebrel_builder.h"

#include "ebrel_setup.h"
#include "species_plan.h"
#include "optimisation_utils.h"  // compute_distance_weights (or wherever it lives)

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>   // std::move

// ------------------------- Local helper functions ----------------------------

namespace {

inline void ensure(bool ok, const std::string& msg) {
  if (!ok) {
    throw std::invalid_argument(msg);
  }
}

} // namespace

//---------------------------- Main function -----------------------------------


BuiltEBREL build_ebrel(
    std::vector<int8_t>  E,
    std::vector<double>  C,
    std::vector<double>  SD,
    std::vector<int>     D,
    std::vector<double>  SxH,
    std::vector<double>  O,
    std::vector<float>   action_weight,
    std::vector<uint8_t> LM,
    std::vector<int8_t>  X0,
    int dim_x,
    int dim_y,
    int n_actions,
    int n_species,
    double sentinel,
    double sigma_in,
    int universal_disp_thres,
    int max_disp_steps,
    int roi_cap,
    bool precompute_W
) {
  
  // ---- BuildEBREL object ----
  BuiltEBREL b;
  
  // ---- Define number of habitat creation options (last action is improvement) ----
  int n_habitats = n_actions - 1;
  
  // ---- Checks on dimensions ----
  ensure(dim_x > 0 && dim_y > 0, "dim_x and dim_y must be positive");
  ensure(n_actions > 0 && n_species > 0, "n_actions and n_species must be positive");
  
  const std::size_t cells = static_cast<std::size_t>(dim_x) * static_cast<std::size_t>(dim_y);
  const std::size_t Asz   = cells * static_cast<std::size_t>(n_actions);
  const std::size_t Ssz   = cells * static_cast<std::size_t>(n_species);
  const std::size_t SxHsz = static_cast<std::size_t>(n_habitats) * static_cast<std::size_t>(n_species);
  
  ensure(E.size()  == cells, "E size mismatch (expected dim_x*dim_y)");
  ensure(LM.size() == cells, "LM size mismatch (expected dim_x*dim_y)");
  ensure(C.size()  == Asz,   "C size mismatch (expected dim_x*dim_y*n_actions)");
  ensure(SD.size() == Ssz,   "SD size mismatch (expected dim_x*dim_y*n_species)");
  ensure(SxH.size()== SxHsz, "SxH size mismatch (expected n_habitats*n_species)");
  ensure(O.size()  == static_cast<std::size_t>(n_species), "O size mismatch (expected n_species)");
  ensure(D.size()  == static_cast<std::size_t>(n_species), "D size mismatch (expected n_species)");
  ensure(action_weight.size() == static_cast<std::size_t>(n_actions), "Action_weight must have length n_actions");
  
  // ---- range checks on E ----
  for (const int8_t v : E) {
    ensure(v >= -1 && v < n_habitats, "E contains value outside [-1, n_habitats-1]");
  }
  
  // ---- Checks on X0 if supplied ----
  const bool has_X0 = !X0.empty();
  
  if (has_X0) {
    ensure(
      X0.size() == cells,
      "X0 size mismatch (expected dim_x*dim_y)"
    );
    
    for (const int8_t action : X0) {
      ensure(
        action >= -1 && action < n_actions,
        "X0 contains value outside [-1, n_actions-1]"
      );
    }
  }
  
  // ---- fill core fields ----
  b.in.dim_x = dim_x;
  b.in.dim_y = dim_y;
  b.in.n_actions = n_actions;
  b.in.n_species = n_species;
  b.in.n_habitats = n_habitats;
  
  b.in.universal_disp_thres = universal_disp_thres;
  b.in.max_disp_steps       = max_disp_steps;
  b.in.roi_cap              = roi_cap;
  
  b.in.E  = std::move(E);
  b.in.C  = std::move(C);
  b.in.SD = std::move(SD);
  b.in.D  = std::move(D);
  b.in.SxH= std::move(SxH);
  b.in.O  = std::move(O);
  b.in.action_weight = std::move(action_weight);
  b.in.LM = std::move(LM);
  b.in.X0 = std::move(X0);
  
  // ---- setup derived state (U, spans, coords, tiles, rowruns, sigma) ----
  create_ebrel_class_object(
    b.in.E, b.in.C, b.in.SD, b.in.D, b.in.SxH, b.in.O, b.in.action_weight, b.in.LM,
    b.in.universal_disp_thres,
    dim_x, dim_y, 
    n_actions, n_species, n_habitats,
    sentinel,
    sigma_in,
    b.sigma,
    b.in.U,
    b.in.row_first_land, b.in.row_last_land,
    b.in.col_first_land, b.in.col_last_land,
    b.in.cell_r, b.in.cell_c,
    b.in.Etiles_per_h,
    b.in.rowruns_cache
  );
  
  // ---- validate supplied X0 against action availability ----
  if (has_X0) {
    for (std::size_t cell = 0; cell < cells; ++cell) {
      const int action = static_cast<int>(b.in.X0[cell]);
      
      // -1 means no action.
      if (action < 0) continue;
      
      const std::size_t idx =
        static_cast<std::size_t>(action) * cells + cell;
      
      ensure(
        b.in.U[idx] == 0u,
        "X0 selects an unavailable action at cell " +
          std::to_string(cell) +
          ", action " +
          std::to_string(action)
      );
    }
  }
  
  // ---- species plan ----
  b.in.species_plan = build_species_plan(
    b.in.SD, b.in.SxH, b.in.D, b.in.O, b.in.action_weight,
    n_actions, n_species, n_habitats, 
    dim_x, dim_y,
    b.in.universal_disp_thres,
    b.in.max_disp_steps,
    b.in.roi_cap,
    b.in.row_first_land, b.in.row_last_land,
    b.in.col_first_land, b.in.col_last_land,
    b.in.E,
    b.in.cell_r, b.in.cell_c
  );
  
  // ---- Initialise improve_count ----
  if (has_X0) {
    b.in.improve_count = initialise_improve_count(
      b.in.X0,
      b.in.species_plan
    );
  } else {
    b.in.improve_count.assign(
      static_cast<std::size_t>(n_species),
      0
    );
  }
  
  // ---- proposal weights W (distance decay) ----
  if (precompute_W) {
    b.in.W = compute_distance_weights(b.in.E, b.in.U, b.in.n_habitats, dim_x, dim_y, b.sigma);
  } else {
    b.in.W.clear();
  }
  
  return b;
}

