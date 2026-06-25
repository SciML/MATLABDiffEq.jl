using SciMLTesting, MATLABDiffEq, JET, Test

run_qa(
    MATLABDiffEq;
    explicit_imports = true,
    ei_kwargs = (;
        # SciMLBase-owned names accessed through their canonical re-exporter
        # DiffEqBase, and the Symbolics-owned MATLABTarget accessed through
        # ModelingToolkit; non-public in the re-exporter, so ignore until they
        # are declared public upstream.
        all_qualified_accesses_via_owners = (;
            ignore = (
                :AbstractODEAlgorithm, :AbstractODEProblem, :__solve,
                :build_solution, :MATLABTarget,
            ),
        ),
        all_qualified_accesses_are_public = (;
            ignore = (
                :AbstractODEAlgorithm, :AbstractODEProblem, :__solve,
                :build_solution, :Stats, :MATLABTarget,
            ),
        ),
    ),
    # no_implicit_imports: the module deliberately `@reexport using DiffEqBase`
    # and `using MATLAB`/`ModelingToolkit`/`PrecompileTools`; making every name
    # explicit is a large, risky refactor against heavy deps. Tracked in
    # https://github.com/SciML/MATLABDiffEq.jl/issues/85
    ei_broken = (:no_implicit_imports,),
)
