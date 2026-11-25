def _defer_runtime_import():
    import sys
    from importlib.util import LazyLoader, find_spec, module_from_spec

    # Register the runtime before exporting the same-named Python API callable.
    runtime_name = f"{__name__}.clipkit"
    if runtime_name in sys.modules:
        return
    spec = find_spec(runtime_name)
    runtime = module_from_spec(spec)
    sys.modules[spec.name] = runtime
    LazyLoader(spec.loader).exec_module(runtime)


_defer_runtime_import()
__all__ = ["clipkit"]


def __getattr__(name):
    if name == "clipkit":
        from .api import clipkit

        globals()[name] = clipkit
        return clipkit
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | set(__all__))
