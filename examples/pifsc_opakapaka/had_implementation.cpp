// Standalone example definition for Quadra's header-only HAD graph pointer.
//
// core/had_quadra.hpp declares:
//   extern threadDefine ADGraph *g_ADGraph;
//
// This file provides the one translation-unit definition needed when compiling
// this example directly with c++ rather than through a larger build target.
#include "../../core/had_quadra.hpp"

namespace had {
threadDefine ADGraph *g_ADGraph = 0;
} // namespace had
