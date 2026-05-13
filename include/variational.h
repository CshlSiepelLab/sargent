#pragma once

#include "sargent.h"

/* Infer fixed effects (beta) and h2_fixed via variational coordinate descent.
   Allocates sg->beta if not already set. */
void var_fixed(sargent_t *sg);
