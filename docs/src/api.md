# API

The MATLABDiffEq algorithm markers are documented under
[Public API](index.md#Public-API) on the home page.

## Reexported SciML common interface

`using MATLABDiffEq` also brings in the parts of the SciML common interface needed to
build an ODE problem, solve it, and inspect the result, so they do not have to be
imported separately. MATLABDiffEq does not define these names -- they are owned and
documented by [SciMLBase](https://docs.sciml.ai/SciMLBase/stable/), and that is where
their documentation lives:

  - Problems: `ODEProblem`, `EnsembleProblem`
  - Functions: `ODEFunction`
  - Solutions: `ODESolution`, `EnsembleSolution`, `EnsembleSummary`, `DEStats`
  - Ensemble algorithms: `EnsembleSerial`, `EnsembleThreads`, `EnsembleDistributed`,
    `EnsembleSplitThreads`, and the `EnsembleAnalysis` module
  - Solving: `solve`, `remake`
  - Return status: `ReturnCode`, `successful_retcode`
  - `NullParameters`

Note that the MATLAB algorithms themselves are public but *not* exported, so they are
still written qualified: `MATLABDiffEq.ode45()`.

Anything else from SciMLBase must be imported from SciMLBase directly. Three groups are
deliberately absent:

  - **DAE, SDE, DDE and every other non-ODE problem type.** `SciMLBase.__solve` is
    defined here only for `AbstractODEProblem`. (`ode15i` names MATLAB's implicit solver,
    but it is still reached through an `ODEProblem`.)
  - **Callbacks.** `__solve` errors on a `callback` keyword ("Callbacks are not supported
    in MATLABDiffEq.jl"), so `ContinuousCallback` and friends are not part of its surface.
  - **The integrator interface** (`init`, `step!`, `solve!`, `reinit!`, ...).
    MATLABDiffEq implements `SciMLBase.__solve` only; it has no integrator.
