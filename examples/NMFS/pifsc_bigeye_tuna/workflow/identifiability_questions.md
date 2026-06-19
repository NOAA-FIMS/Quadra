# Bigeye Tuna Identifiability Roadmap

This roadmap describes a diagnostic-guided modeling workflow. The goal is not to
port MULTIFAN-CL. The goal is to build a better model-development process where
complexity is added only when the data and diagnostics support it.

## Core principle

Each model level must answer:

1. What feature was added?
2. What information should identify that feature?
3. What diagnostics improved or degraded?
4. What parameter confounding appeared?
5. Should the feature be retained, simplified, pooled, or removed?

## Level 0: single-region age-structured model

Feature: one region, one aggregate catch series, one index, age structure,
recruitment deviations, and fixed effects for recruitment scale, fishing
mortality, and catchability.

Question: Can recruitment scale, fishing mortality, catchability, and
recruitment deviations be separated with one catch series and one abundance
index?

Diagnostics:

- convergence status
- fixed-effect gradient norm
- Huu positive definiteness
- Huu condition number
- effective sparsity
- effective bandwidth
- recruitment-deviation correlation graph
- parameter influence rankings
- objective component contributions

Decision rule: If q, F, and recruitment scale are strongly confounded, do not
add fleets or movement yet.

## Level 1: multiple fleets

Feature: multiple catch fleets, fleet-specific q or scaling terms, and
fleet-specific selectivity where justified.

Question: Which fleets identify abundance, fishing mortality, catchability, and
selectivity?

Decision rule: Fleet-specific q/selectivity should be retained only when
diagnostics show that the fleet contributes distinct information.

## Level 2: fleet-specific composition data

Feature: length or age composition by fleet.

Question: Which composition data identify selectivity, and which are redundant
or conflicting?

Decision rule: If a fleet's composition data create high curvature but poor fit
or strong confounding, consider pooling selectivity or simplifying bins.

## Level 3: spatial structure

Feature: abundance by region and region-specific index/fleet mapping.

Question: Are regional abundance trends identified separately?

Decision rule: Do not add movement until regional abundance is at least partly
identifiable.

## Level 4: movement

Feature: movement among regions.

Question: Can movement be distinguished from regional recruitment, regional
abundance, and catchability?

Decision rule: Movement parameters should be simplified, pooled, or regularized
if the data cannot identify them.

## Level 5: tagging

Feature: releases, recaptures, reporting rates, and movement likelihood.

Question: Does tagging identify movement, reporting rate, or both?

Decision rule: Retain tagging complexity only if it reduces movement uncertainty
or clearly separates movement from reporting.
