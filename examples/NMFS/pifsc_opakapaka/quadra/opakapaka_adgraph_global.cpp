// Defines the global/thread-local HAD graph pointer for standalone Quadra
// examples.
//
// Several test binaries already link an implementation translation unit.
// The opakapaka fair benchmark is built directly with c++, so it needs one too.
#include "../../../../core/had_quadra.hpp"

namespace had {
threadDefine ADGraph *g_ADGraph = nullptr;
} // namespace had
