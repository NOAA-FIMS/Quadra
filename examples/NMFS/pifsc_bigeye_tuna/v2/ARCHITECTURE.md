# Composable Assessment Architecture (CAA)

Matthew Supernaw

NOAA Fisheries

Version 0.1


*"Small Steps. Explicit State. Composable Science."*

## Core Principles

State owns memory.

Steps own algorithms.

## The Eight Rules

1. State owns memory.

2. Steps own algorithms.

3. Data is immutable.

4. Parameters describe the model.

5. AssessmentCycle owns execution order.

6. Every scientific concept has one owner.

7. Every Step is independently testable.

8. Composition replaces inheritance.

Every contribution to CAA should be evaluated against these principles before it is evaluated against the implementation.

------------------------------------------------------------

## Overview

The Composable Assessment Architecture (CAA) is a reference architecture for
building stock assessment models from small, reusable scientific components.

Rather than implementing an entire assessment model as one large class,
an assessment is built by composing independent scientific steps.

Each step performs one scientific calculation.

The framework is responsible for orchestration.

Scientists are responsible for equations.

---

# Philosophy

An assessment model consists of five concepts.

```
Assessment
    ├── Data
    ├── Parameters
    ├── State
    ├── Steps
    └── Assessment Cycle
```

These concepts are intentionally independent.

---

# Data

Data represents information provided to the assessment.

Examples include

- Catch observations
- Survey observations
- Age compositions
- Length compositions
- Environmental covariates
- Fleet definitions
- Population definitions

Data is immutable.

Scientific steps never modify data.

---

# Parameters

Parameters describe the model.

Parameters may be

- estimated
- fixed
- shared
- hierarchical

Examples include

- recruitment parameters
- selectivity parameters
- fishing mortality parameters
- movement parameters

Parameters describe the model.

They do not perform calculations.

---

# State

State contains quantities produced while executing the assessment.

Examples include

- weight-at-age
- maturity-at-age
- natural mortality-at-age
- numbers-at-age
- spawning biomass
- fishing mortality-at-age
- catch-at-age
- predicted observations
- likelihood components

State is explicit.

There is no hidden state.

---

# Scientific Steps

A scientific step performs exactly one scientific calculation.

Examples include

Life History

- Logistic Maturity
- Lorenzen Mortality
- Weight-at-Age

Population

- Recruitment
- Survival
- Plus Group

Movement

- Annual Movement
- Seasonal Movement

Fleet

- Logistic Selectivity
- Fishing Mortality
- Baranov Catch

Observation

- Biomass Index
- Age Composition

Likelihood

- Catch Likelihood
- Survey Likelihood
- Age Composition Likelihood

Every step

- has one responsibility
- has explicit inputs
- has explicit outputs
- is independently executable
- is independently testable

---

# Assessment Cycle

The assessment cycle orchestrates scientific steps.

For a typical annual assessment,

```
Life History

↓

for each population

    Population Dynamics

↓

Movement

↓

for each fleet

    Selectivity

    Fishing Mortality

    Catch

↓

Observation Models

↓

Likelihood

↓

Objective Function

↓

Optimization
```

The assessment cycle owns execution order.

Scientific steps own scientific calculations.

---

# Design Principles

## Composition over inheritance

Models are assembled from reusable scientific steps.

Framework behavior should emerge from composition,
not deep inheritance hierarchies.

---

## Explicit state

State should always be visible.

Scientific steps receive the state they require and
produce the state they own.

---

## Infrastructure orchestrates

The framework owns

- execution order
- optimization
- automatic differentiation
- Laplace approximation
- scheduling
- parallelization

The scientific model owns

- biological equations
- observation equations
- likelihood equations

---

## One scientific concept per component

Every component should answer exactly one scientific question.

Examples

```
LogisticSelectivity

BaranovCatch

LorenzenMortality

AnnualMovement
```

Avoid components that perform multiple unrelated calculations.

---

## Regression-first development

Every scientific step should have an executable regression test.

Complex models are built by composing independently verified components.

---

# Reference Architecture

```
architecture/

    data/

    parameters/

    state/

    steps/

        life_history/

        population/

        movement/

        fleet/

        observation/

        likelihood/

    assessment/

        annual_cycle.hpp

        objective.hpp

        optimizer.hpp
```

---

# Goals

The Composable Assessment Architecture is designed to

- simplify assessment development
- improve readability
- encourage code reuse
- reduce coupling
- simplify testing
- separate infrastructure from science
- make assessment models easier to teach
- support multiple optimization backends
- support multiple automatic differentiation libraries

---

# Vision

Assessment models should read like assessment science.

A scientist should be able to understand the structure of a model
without understanding the implementation details.

The code should explain the science.

The framework should disappear.

The purpose of CAA is to let scientists focus on science while the framework handles infrastructure.
