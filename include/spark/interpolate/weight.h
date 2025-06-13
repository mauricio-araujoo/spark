#pragma once

#include "spark/particle/species.h"
#include "spark/threads/pool.h"

namespace spark::interpolate {

template <class GridType, unsigned NX, unsigned NV>
void weight_to_grid(const spark::particle::ChargedSpecies<NX, NV>& species, GridType& out, spark::threads::ThPool& pool);

}
