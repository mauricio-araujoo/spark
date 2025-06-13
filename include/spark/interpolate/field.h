#pragma once

#include "spark/threads/pool.h"
#include "spark/particle/species.h"
#include "spark/spatial/grid.h"

namespace spark::interpolate {

template <typename T, unsigned NX, unsigned NV>
void field_at_particles(const spatial::TUniformGrid<T, NX>& field,
                        const particle::ChargedSpecies<NX, NV>& species,
                        core::TMatrix<T, 1>& out,
                        spark::threads::ThPool& pool);
}