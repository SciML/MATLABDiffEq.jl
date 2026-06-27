using SciMLTesting, MATLABDiffEq, JET, Test

run_qa(
    MATLABDiffEq;
    explicit_imports = true,
    ei_kwargs = (;
        # DiffEqBase.__solve (SciMLBase-owned) is the documented solver extension
        # point re-exported by DiffEqBase, and the Symbolics-owned MATLABTarget is
        # reached through ModelingToolkit; both are accessed through a re-exporter
        # that is not the owner.
        all_qualified_accesses_via_owners = (;
            ignore = (
                :__solve, :MATLABTarget,
            ),
        ),
        # Still non-public upstream: __solve (SciMLBase) and MATLABTarget
        # (Symbolics). Drop once they are declared public.
        all_qualified_accesses_are_public = (;
            ignore = (
                :__solve, :MATLABTarget,
            ),
        ),
    ),
    # no_implicit_imports: the module deliberately `@reexport using DiffEqBase`
    # and `using MATLAB`/`ModelingToolkit`/`PrecompileTools`; making every name
    # explicit is a large, risky refactor against heavy deps. Tracked in
    # https://github.com/SciML/MATLABDiffEq.jl/issues/85
    ei_broken = (:no_implicit_imports,),
)
