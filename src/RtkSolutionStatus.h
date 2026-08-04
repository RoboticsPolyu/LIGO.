#pragma once

#include <cstddef>
#include <cstdint>

enum class RtkSolutionType : uint8_t
{
  NoAmbiguity = 0,
  Float = 1,
  Fixed = 2
};

// A navigation state is considered fixed whenever at least one active
// ambiguity is constrained to an integer. This takes precedence over any
// remaining float ambiguities because the exported state already depends on
// an integer constraint.
inline RtkSolutionType classifyRtkSolution(size_t active_float,
                                           size_t active_fixed)
{
  if (active_fixed > 0) return RtkSolutionType::Fixed;
  if (active_float > 0) return RtkSolutionType::Float;
  return RtkSolutionType::NoAmbiguity;
}

inline const char *rtkSolutionTypeName(RtkSolutionType type)
{
  switch (type)
  {
    case RtkSolutionType::NoAmbiguity: return "NO_AMBIGUITY";
    case RtkSolutionType::Float: return "FLOAT";
    case RtkSolutionType::Fixed: return "FIXED";
  }
  return "NO_AMBIGUITY";
}
